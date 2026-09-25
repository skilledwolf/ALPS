/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2011 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#ifndef PARAPACK_EXCHANGE_MULTI_H
#define PARAPACK_EXCHANGE_MULTI_H

#include "exchange.h"

namespace alps {
namespace parapack {
namespace exmc {

template<typename WALKER, typename INITIALIZER> struct multiple_initializer_helper;

template<typename WALKER>
struct multiple_initializer_helper<WALKER, no_initializer> {
  typedef WALKER walker_type;
  static boost::shared_ptr<walker_type>
  create_walker(boost::mpi::communicator const& comm, alps::Parameters const& params, no_initializer const&) {
    return boost::shared_ptr<walker_type>(new walker_type(comm, params));
  }

  static void init_observables(boost::shared_ptr<walker_type> walker_ptr,
    alps::Parameters const& params, no_initializer const&, alps::ObservableSet& obs) {
    walker_ptr->init_observables(params, obs);
  }
  static void run_walker(boost::shared_ptr<walker_type> walker_ptr, no_initializer const&, alps::ObservableSet& obs) {
    walker_ptr->run(obs);
  }
};

template<typename WALKER, typename INITIALIZER>
struct multiple_initializer_helper {
  typedef WALKER walker_type;
  typedef INITIALIZER initializer_type;
  static boost::shared_ptr<walker_type>
  create_walker(boost::mpi::communicator const& comm, alps::Parameters const& params, initializer_type const& init) {
    return boost::shared_ptr<walker_type>(new walker_type(comm, params, init));
  }
  static void init_observables(boost::shared_ptr<walker_type> walker_ptr, alps::Parameters const& params, initializer_type const& init, alps::ObservableSet& obs) {
    walker_ptr->init_observables(params, init, obs);
  }
  static void run_walker(boost::shared_ptr<walker_type> walker_ptr, initializer_type const& init, alps::ObservableSet& obs) {
    walker_ptr->run(init, obs);
  }
};

} // end namespace exmc
  
template<typename WALKER, typename INITIALIZER = exmc::no_initializer>
class multiple_parallel_exchange_worker : public mc_worker {
private:
  friend struct exmc::exchange_algorithm;
  typedef mc_worker super_type;
  typedef WALKER walker_type;
  typedef typename walker_type::weight_parameter_type weight_parameter_type;
  typedef INITIALIZER initializer_type;
  typedef exmc::multiple_initializer_helper<walker_type, initializer_type> helper;
  typedef exmc::walker_direc walker_direc;

public:
  static std::string version() { return walker_type::version(); }
  static void print_copyright(std::ostream& out) { walker_type::print_copyright(out); }

  multiple_parallel_exchange_worker(boost::mpi::communicator const& comm, alps::Parameters const& params)
    : super_type(params), process_helper_(comm, params.value_or_default("PROCESS_PER_WORKER", 1)),
      comm_(process_helper_.comm_ctrl()), work_(process_helper_.comm_work()),
      head_(process_helper_.comm_head()), init_(params), beta_(params), mcs_(params),
      num_returnee_(0) {

    if (head_ && head_.rank() == 0) {
      if (comm_.rank() != 0) {
        std::cerr << "Error: total number of processes should be a multiple of PROCESS_PER_WORKER\n";
        boost::throw_exception(std::runtime_error(
          "total number of processes should be a multiple of PROCESS_PER_WORKER"));
      }
    }
    
    if (head_) {
      boost::tie(nrep_local_, offset_local_) = calc_nrep(head_.rank());
      if (nrep_local_ == 0) {
        std::cerr << "Error: number of replicas is smaller than number of processes\n";
        boost::throw_exception(std::runtime_error(
          "number of replicas is smaller than number of processes"));
      }
    }
    broadcast(work_, nrep_local_, 0);
    broadcast(work_, offset_local_, 0);
    if (head_ && head_.rank() == 0) exmc::exchange_algorithm::init_partition(*this, head_.size());
    exmc::exchange_algorithm::init_walkers(*this, params, tid_local_, nrep_local_, offset_local_,
      [&](alps::Parameters const& wp) { return helper::create_walker(work_, wp, init_); });
    if (head_ && head_.rank() == 0) exmc::exchange_algorithm::init_state(*this);
    if (mcs_.exchange() && head_) wp_local_.resize(nrep_local_);
  }
  virtual ~multiple_parallel_exchange_worker() {}

  void init_observables(alps::Parameters const& params, std::vector<alps::ObservableSet>& obs) {
    int nrep = beta_.size();
    obs.resize(nrep);
    for (int p = 0; p < nrep; ++p)
      helper::init_observables(walker_[0], params, init_, obs[p]);
    if (head_ && head_.rank() == 0) {
      for (int p = 0; p < nrep; ++p) {
        exmc::exchange_algorithm::init_observables(obs[p], p, nrep, mcs_.exchange());
      }
      if (mcs_.exchange()) obs[0] << SimpleRealObservable("EXMC: Average Inverse Round-Trip Time");
    }
  }

  void run(std::vector<alps::ObservableSet>& obs) {
    ++mcs_;

    int nrep = beta_.size();

    if (head_ && head_.rank() == 0) {
      for (int p = 0; p < nrep; ++p) {
        add_constant(obs[p]["EXMC: Temperature"], 1. / beta_[p]);
        add_constant(obs[p]["EXMC: Inverse Temperature"], beta_[p]);
      }
    }

    // MC update of each replica
    for (int w = 0; w < nrep_local_; ++w) {
      int p = tid_local_[w];
      walker_[w]->set_beta(beta_[p]);
      helper::run_walker(walker_[w], init_, obs[p]);
    }

    // replica exchange process
    if (mcs_.exchange() && (mcs_() % mcs_.interval()) == 0) {

      bool continue_stage = false;
      bool next_stage = false;

      if (head_) {
        for (int w = 0; w < nrep_local_; ++w) wp_local_[w] = walker_[w]->weight_parameter();
        if (head_.rank() == 0) {
          std::copy(wp_local_.begin(), wp_local_.begin() + nrep_local_, wp_.begin());
          for (int p = 1; p < nreps_.size(); ++p)
            head_.recv(p, 0, &wp_[offsets_[p]], nreps_[p]);
          for (int w = 0; w < nrep; ++w) {
            int p = tid_[w];
            weight_parameters_[p] += wp_[w];
          }
        } else {
          head_.send(0, 0, &wp_local_[0], nrep_local_);
        }

        if (head_.rank() == 0) {
          exmc::exchange_algorithm::exchange<true>(*this, obs);

          if (mcs_.doing_optimization() && mcs_.stage_count() == mcs_.stage_sweeps()) {

            exmc::exchange_algorithm::optimize(*this, obs,
              [&] { next_stage = true; }, [&] { continue_stage = true; }, 5);

            // check whether all the replicas have revisited the highest temperature or not
            exmc::exchange_algorithm::check_thermalization(*this, [&] { continue_stage = true; });
          }
        }
      
        // broadcast EXMC results
        broadcast(head_, continue_stage, 0);
        broadcast(head_, next_stage, 0);
      }
      broadcast(work_, continue_stage, 0);
      broadcast(work_, next_stage, 0);
      if (continue_stage) mcs_.continue_stage();
      if (next_stage) mcs_.next_stage();
      if (head_) {
        if (head_.rank() == 0) {
          for (int p = 1; p < nreps_.size(); ++p)
            head_.send(p, 0, &tid_[offsets_[p]], nreps_[p]);
          std::copy(tid_.begin(), tid_.begin() + nrep_local_, tid_local_.begin());
        } else {
          head_.recv(0, 0, &tid_local_[0], nrep_local_);
        }
      }
      for (int i = 0; i < nrep_local_; ++i)
        broadcast(work_, tid_local_[i], 0);
    }
  }

  void save(alps::ODump& dp) const {
    dp << beta_ << mcs_ << tid_local_;
    if (comm_.rank() == 0) dp << tid_ << wid_ << direc_ << num_returnee_ << weight_parameters_;
    for (int i = 0; i < nrep_local_; ++i) walker_[i]->save(dp);
  }
  void load(alps::IDump& dp) {
    dp >> beta_ >> mcs_ >> tid_local_;
    if (comm_.rank() == 0) dp >> tid_ >> wid_ >> direc_ >> num_returnee_ >> weight_parameters_;
    for (int i = 0; i < nrep_local_; ++i) walker_[i]->load(dp);
  }

  bool is_thermalized() const { return mcs_.is_thermalized(); }
  double progress() const { return mcs_.progress(); }

  static void evaluate_observable(alps::ObservableSet& obs) {
    walker_type::evaluate_observable(obs);
  }

protected:
  std::pair<int, int> calc_nrep(int id) const {
    // Replicas are distributed across worker groups, not individual ranks.
    return exmc::exchange_algorithm::partition(beta_.size(), head_.size(), id);
  }

private:
  alps::process_helper_mpi process_helper_;
  boost::mpi::communicator comm_, work_, head_;

  int nrep_local_;           // number of walkers (replicas) on this process
  int offset_local_;         // first (global) id of walker on this process
  int nrep_max_;             // [master only] maximum number of walkers (replicas) on a process
  std::vector<int> nreps_;   // [master only] number of walkers (replicas) on each process
  std::vector<int> offsets_; // [master only] first (global) id of walker on each process

  initializer_type init_;
  std::vector<boost::shared_ptr<walker_type> > walker_; // [0..nrep_local_)

  exmc::inverse_temperature_set beta_;
  exmc::exchange_steps mcs_;
  std::vector<int> tid_local_; // temperature id of each walker (replica)
  std::vector<int> tid_;       // [master only] temperature id of each walker (replica)
  std::vector<int> wid_;       // [master only] walker (replica) id at each temperature
  std::vector<int> direc_;     // [master only] direction of each walker (replica)
  int num_returnee_;           // [master only] number of walkers returned to highest temperature
  std::vector<weight_parameter_type> weight_parameters_; // [master only]

  // working space
  std::vector<weight_parameter_type> wp_local_;
  std::vector<weight_parameter_type> wp_; // [master only]
  std::vector<double> upward_;            // [master only]
  std::vector<double> accept_;            // [master only]
  std::vector<int> permutation_;          // [master only]
};

} // end namespace parapack
} // end namespace alps

#endif // PARAPACK_EXCHANGE_MULTI_H
