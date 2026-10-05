/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2010 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#ifndef PARAPACK_EXAMPLE_MULTIPLE_ISING_H
#define PARAPACK_EXAMPLE_MULTIPLE_ISING_H

#include <alps/parapack/worker.h>
#include <cmath>
#include <functional>
#include <limits>

namespace mpi = boost::mpi;

// a vector with 'tabs' on both sides: elements at i = -1 as well as
// i = size() are also accessible
template<typename T>
class tabbed_vector {
public:
  typedef T value_type;
  typedef typename std::vector<value_type>::size_type size_type;

  explicit tabbed_vector(size_type n = 0) : vector_(n + 2) {}
  explicit tabbed_vector(size_type n, value_type x) : vector_(n + 2, x) {}
  void resize(size_type n, value_type x = value_type()) { vector_.resize(n + 2, x); }

  size_type size() const { return vector_.size() - 2; }
  value_type const& operator[](int i) const { return vector_[i+1]; }
  value_type& operator[](int i) { return vector_[i+1]; }

  void save(alps::ODump& dp) const { dp << vector_; }
  void load(alps::IDump& dp) { dp >> vector_; }

private:
  std::vector<value_type> vector_;
};

template<typename T>
alps::ODump& operator<<(alps::ODump& dp, tabbed_vector<T> const& v) { v.save(dp); return dp; }

template<typename T>
alps::IDump& operator>>(alps::IDump& dp, tabbed_vector<T>& v) { v.load(dp); return dp; }


class parallel_ising_worker : public alps::parapack::mc_worker {
private:
  typedef alps::parapack::mc_worker super_type;

public:
  parallel_ising_worker(mpi::communicator const& comm, alps::Parameters const& params)
    : super_type(params), comm_(comm), beta_(std::numeric_limits<double>::quiet_NaN()),
      coupling_(1), length_(0), loclen_(0), mcs_(params), energy_(0) {
    // Exchange workers assign beta after construction. Invalid local input must
    // be rejected on every rank before any rank enters a halo exchange.
    bool valid = true;
    try {
      if (params.defined("T")) {
        const double temperature = evaluate("T", params);
        beta_ = 1 / temperature;
        valid = temperature > 0 && std::isfinite(temperature) && std::isfinite(beta_);
      }
      coupling_ = params.defined("J") ? evaluate("J", params) : 1.0;
      const double length = evaluate("L", params);
      valid = valid && std::isfinite(coupling_) && std::isfinite(length) &&
        length >= 2.0 * comm_.size() && length <= std::numeric_limits<int>::max() &&
        length == std::floor(length);
      if (valid) length_ = static_cast<int>(length);
    } catch (std::exception const&) {
      valid = false;
    }
    double model[]{double(length_), coupling_, std::isnan(beta_) ? -1 : beta_};
    mpi::broadcast(comm_, model, 3, 0);
    valid = valid && model[0] == length_ && model[1] == coupling_ &&
      model[2] == (std::isnan(beta_) ? -1 : beta_);
    if (!mpi::all_reduce(comm_, valid, std::logical_and<bool>()))
      throw std::invalid_argument("Spatial Ising requires matching model inputs, finite J, positive finite T when specified, and integer L >= 2 * ranks");
    // system size and local system size
    if (comm_.rank() == 0)
      loclen_ = length_ - (comm_.size() - 1) * (length_ / comm_.size());
    else
      loclen_ = length_ / comm_.size();

    // configuration
    spins_.resize(loclen_);
    for (int i = 0; i < loclen_; ++i) spins_[i] = (uniform_01() < 0.5 ? 1 : 0);
    copy2right();
    copy2left();
  }
  virtual ~parallel_ising_worker() {}

  void init_observables(alps::Parameters const&, alps::ObservableSet& obs) {
    if (comm_.rank() == 0)
      obs << alps::SimpleRealObservable("Temperature")
          << alps::SimpleRealObservable("Inverse Temperature")
          << alps::SimpleRealObservable("Number of Sites")
          << alps::RealObservable("Energy")
          << alps::RealObservable("Energy^2")
          << alps::RealObservable("Magnetization")
          << alps::RealObservable("Magnetization^2")
          << alps::RealObservable("Magnetization^4");
  }

  bool is_thermalized() const { return mcs_.is_thermalized(); }
  double progress() const { return mcs_.progress(); }

  void run(alps::ObservableSet& obs) {
    const double minimum_beta = mpi::all_reduce(comm_, std::isfinite(beta_) && beta_ >= 0 ? beta_ : -1., mpi::minimum<double>());
    if (minimum_beta < 0)
      throw std::invalid_argument("Set a finite nonnegative inverse temperature before updating spatial Ising");
    if (minimum_beta != mpi::all_reduce(comm_, beta_, mpi::maximum<double>()))
      throw std::invalid_argument("Spatial Ising requires the same inverse temperature on every rank");
    ++mcs_;

    for (int i = 0; i < loclen_; ++i) {
      // H = -J sum(s_i s_{i+1}), with s_i = 2 * spins_[i] - 1.
      double diff = coupling_ * (4 - 4 * ((spins_[i-1] ^ spins_[i]) + (spins_[i] ^ spins_[i+1])));
      if (uniform_01() < 0.5 * (1 + std::tanh(-0.5 * beta_ * diff))) spins_[i] ^= 1;
      if (i == 0) copy2left();
      if (i == loclen_ - 1) copy2right();
    }

    // measurements
    energy_ = 0;
    double mag = 0;
    for (int i = 0; i < loclen_; ++i) {
      energy_ += coupling_ * (2 * (spins_[i] ^ spins_[i+1]) - 1);
      mag += (2 * spins_[i] - 1);
    }
    if (comm_.rank() == 0) {
      double energy_out, mag_out;
      reduce(comm_, energy_, energy_out, std::plus<double>(), 0);
      reduce(comm_, mag, mag_out, std::plus<double>(), 0);
      energy_ = energy_out;
      mag = mag_out;
    } else {
      reduce(comm_, energy_, std::plus<double>(), 0);
      reduce(comm_, mag, std::plus<double>(), 0);
    }

    if (comm_.rank() == 0) {
      add_constant(obs["Temperature"], 1/beta_);
      add_constant(obs["Inverse Temperature"], beta_);
      add_constant(obs["Number of Sites"], (double)length_);
      obs["Energy"] << energy_;
      obs["Energy^2"] << energy_ * energy_;
      obs["Magnetization"] << mag;
      obs["Magnetization^2"] << mag * mag;
      obs["Magnetization^4"] << mag * mag * mag * mag;
    }
  }

  // for exchange Monte Carlo
  typedef double weight_parameter_type;
  void set_beta(double beta) { beta_ = beta; } // Validation is collective at the next update.
  weight_parameter_type weight_parameter() const { return energy_; }
  static double log_weight(weight_parameter_type gw, double beta) { return - beta * gw; }

  void save(alps::ODump& dp) const { dp << mcs_ << spins_ << energy_; }
  void load(alps::IDump& dp) { dp >> mcs_ >> spins_ >> energy_; }

protected:
  void copy2right() {
    if (comm_.size() == 1) {
      spins_[-1] = spins_[loclen_-1];
    } else {
      comm_.sendrecv((comm_.rank() + 1) % comm_.size(), 0, spins_[loclen_-1],
                    (comm_.rank() + comm_.size() - 1) % comm_.size(), 0, spins_[-1]);
    }
  }

  void copy2left() {
    if (comm_.size() == 1) {
      spins_[loclen_] = spins_[0];
    } else {
      comm_.sendrecv((comm_.rank() + comm_.size() - 1) % comm_.size(), 0, spins_[0],
                    (comm_.rank() + 1) % comm_.size(), 0, spins_[loclen_]);
    }
  }

private:
  mpi::communicator comm_;

  // parameteters
  double beta_;
  double coupling_;
  int length_;
  int loclen_;

  // configuration (need checkpointing)
  alps::mc_steps mcs_;
  tabbed_vector<int> spins_;
  double energy_;
};

#endif // PARAPACK_EXAMPLE_MULTIPLE_ISING_H
