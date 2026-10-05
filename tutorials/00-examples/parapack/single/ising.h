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

#pragma once
#include <alps/parapack/worker.h>
#include "kernel.hpp"

// Temporary adapter for the remaining legacy exchange and spatial-MPI callers.
class single_ising_worker : public alps::parapack::lattice_mc_worker<> {
    using super_type=alps::parapack::lattice_mc_worker<>;
public:
    single_ising_worker(alps::Parameters const& p) : super_type(p),
        beta_(p.defined("T") ? 1/evaluate("T",p) : 0), mcs_(p),
        state_(static_cast<alps::graph_helper<> const&>(*this),p.defined("J") ? evaluate("J",p) : 1.,[&]{return uniform_01();}) {}
    void init_observables(alps::Parameters const&,alps::ObservableSet& obs) {
        obs<<alps::SimpleRealObservable("Temperature")<<alps::SimpleRealObservable("Inverse Temperature")
           <<alps::SimpleRealObservable("Number of Sites")<<alps::RealObservable("Energy")
           <<alps::RealObservable("Energy^2")<<alps::RealObservable("Magnetization")
           <<alps::RealObservable("Magnetization^2")<<alps::RealObservable("Magnetization^4");
    }
    bool is_thermalized() const {return mcs_.is_thermalized();}
    double progress() const {return mcs_.progress();}
    void run(alps::ObservableSet& obs) {
        ++mcs_; state_.step(beta_,[&]{return uniform_01();});
        const auto x=state_.sample();
        add_constant(obs["Temperature"],1/beta_); add_constant(obs["Inverse Temperature"],beta_);
        add_constant(obs["Number of Sites"],x[0]);
        const char* names[]={"Energy","Energy^2","Magnetization","Magnetization^2","Magnetization^4"};
        for (size_t i=0;i<5;++i) obs[names[i]]<<x[i+1];
    }
    using weight_parameter_type=double;
    void set_beta(double beta) {beta_=beta;}
    double weight_parameter() const {return -state_.energy();}
    static double log_weight(double weight,double beta) {return beta*weight;}
    void save(alps::ODump& out) const {out<<mcs_<<state_.spins()<<state_.energy();}
    void load(alps::IDump& in) {
        std::vector<int> spins; double energy;
        in>>mcs_>>spins>>energy; state_.restore(std::move(spins));
    }
private:
    double beta_;
    alps::mc_steps mcs_;
    ising_kernel state_;
};

class ising_evaluator : public alps::parapack::simple_evaluator {
public:
  ising_evaluator(alps::Parameters const&) {}
  virtual ~ising_evaluator() {}

  void evaluate(alps::ObservableSet& obs) const {
    if (obs.has("Inverse Temperature") && obs.has("Number of Sites") &&
        obs.has("Energy") && obs.has("Energy^2")) {
      alps::RealObsevaluator beta = obs["Inverse Temperature"];
      alps::RealObsevaluator n = obs["Number of Sites"];
      alps::RealObsevaluator ene = obs["Energy"];
      alps::RealObsevaluator ene2 = obs["Energy^2"];
      if (beta.count() && n.count() && ene.count() && ene2.count()) {
        alps::RealObsevaluator c("Specific Heat");
        c = beta.mean() * beta.mean() * (ene2 - ene * ene) / n.mean();
        obs.addObservable(c);
      }
    }
    if (obs.has("Magnetization^2") && obs.has("Magnetization^4")) {
      alps::RealObsevaluator m2 = obs["Magnetization^2"];
      alps::RealObsevaluator m4 = obs["Magnetization^4"];
      if (m2.count() && m4.count()) {
        alps::RealObsevaluator binder("Binder Ratio of Magnetization");
        binder = m2 * m2 / m4;
        obs.addObservable(binder);
      }
    }
    if (obs.has("Inverse Temperature")) obs.erase("Inverse Temperature");
    if (obs.has("Temperature")) obs.erase("Temperature");
  }
};
