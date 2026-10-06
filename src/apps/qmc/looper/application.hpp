// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#pragma once
#include "loop_worker.h"
#include <alps/mc/replica_exchange.hpp>
#include <alps/mc/replica_parallel.hpp>
#ifdef ALPS_HAVE_MPI
#include <boost/serialization/utility.hpp>
#endif
#include <optional>
#include <variant>

namespace looper {
class application {
  using ct_worker=loop_worker<path_integral>;
  using sse_worker=loop_worker<sse>;
  using worker=std::variant<std::unique_ptr<ct_worker>,std::unique_ptr<sse_worker>>;
  using weight=ct_worker::weight_parameter_type;
  struct replica_statistics : native_qmc::simulation {
    replica_statistics(alps::params const& p,size_t bins,size_t chain,native_qmc::simulation const& prototype):simulation(p,bins,chain) {
      labels_=prototype.measurement_labels();signed_names_=prototype.signed_measurements();
      for (auto const& name:alps::mc::batch_names(prototype.get_measurements()))
        alps::mc::add_measurement(*this,name,prototype.measurement(name)->size(),bins);
    }
    void update() override {}
    void measure() override {}
    double fraction_completed() const override {return 0.;}
  };
  // Physical state exists only on its owner. Each temperature retains the
  // same ordered raw observations, independent of where its walker runs.
  std::vector<worker> workers_;
  std::vector<std::unique_ptr<replica_statistics>> stats_;
  std::unique_ptr<alps::mc::replica_checkpoint> checkpoint_;
  alps::params parameters_;
  alps::mc::replica_parallel group_;
  size_t bins_,chain_,sites_;
  std::optional<alps::mc::replica_exchange<weight>> exchange_;

  template<class F> decltype(auto) visit(size_t i,F&& f) const {
    return std::visit([&](auto const& value)->decltype(auto){return f(*value);},workers_.at(i));
  }
  native_qmc::simulation& statistics(size_t i) const {
    if (exchange_) return *stats_.at(i);
    return visit(i,[](auto& value)->native_qmc::simulation& {return value;});
  }
public:
  application(alps::params const& p,size_t bins,size_t chain,alps::mc::replica_parallel group):parameters_(p),group_(std::move(group)),bins_(bins),chain_(chain) {
    auto algorithm=p.value_or<std::string>("ALGORITHM","loop");
    bool ladder=algorithm.find("exchange")!=std::string::npos;
    if (group_.distributed && !ladder) throw std::invalid_argument("parallel = replicas requires an exchange algorithm");
    if (ladder) exchange_.emplace(p,chain,weight::Zero());
    size_t n=ladder ? exchange_->size() : 1;
    workers_.resize(n);
    auto create=[&](size_t i) {
      // Sampling beta is separate from the common XML model parameters.
      double initial_beta=ladder ? exchange_->beta(0) : 0.;
      size_t seed_offset=ladder ? chain*(n+1)+i : chain;
      if (algorithm.find("sse")!=std::string::npos) workers_[i]=std::make_unique<sse_worker>(p,bins,seed_offset,initial_beta);
      else workers_[i]=std::make_unique<ct_worker>(p,bins,seed_offset,initial_beta);
    };
    size_t prototype=n;
    for (size_t i=0;i<n;++i) if (group_.owns_walker(i)) {
      create(i);if (prototype==n) prototype=i;
    }
    // A temporary prototype also validates the model on idle ranks.
    const bool idle=prototype==n;
    if (idle) {prototype=0;create(prototype);}
    sites_=visit(prototype,[](auto const& value){return value.site_count();});
    if (ladder) {
      for (size_t i=0;i<n;++i) stats_.push_back(visit(prototype,[&](auto const& value) {
        return std::make_unique<replica_statistics>(p,bins,chain,value);
      }));
      exchange_->init_diagnostics([&](size_t i,char const* name) {
        statistics(i).add_measurement(name,1,false);
        if (group_.owns_walker(i)) visit(i,[&](auto& value){value.add_measurement(name,1,false);});
      });
    }
    if (idle) workers_[prototype]=std::unique_ptr<ct_worker>{};
  }
  static alps::params checkpoint_parameters(alps::params p) {return native_qmc::simulation::checkpoint_parameters(p);}
  size_t replicas() const {return workers_.size();}
  bool ladder() const {return bool(exchange_);}
  auto temperatures() const {return exchange_->values();}
  void update() {
    if (!exchange_) {visit(0,[](auto& value){value.update();});return;}
    exchange_->step([&](auto const& walkers,auto const& betas,bool sampling) {
      std::vector<native_qmc::simulation::samples_type> samples(replicas());
      group_.update([&] {
        for (size_t i=0;i<replicas();++i) {
          auto& stats=statistics(i);stats.record_measurements(sampling);
          if (group_.owns_walker(walkers[i]))
            samples[i]=stats.sample([&]{visit(walkers[i],[&](auto& value){value.run(stats,betas[i]);});});
        }
      });
      for (size_t i=0;i<replicas();++i) {
        group_.broadcast(samples[i],walkers[i]);
        statistics(i).record(samples[i]);
      }
    },[&] {
      std::vector<double> values(replicas()*2);
      for (size_t w=0;w<replicas();++w) if (group_.owns_walker(w)) {
        auto weight=visit(w,[](auto const& value){return value.weight_parameter();});
        values[2*w]=weight[0];values[2*w+1]=weight[1];
      }
      group_.collect(values);
      std::vector<weight> weights(replicas());
      for (size_t w=0;w<replicas();++w) weights[w]<<values[2*w],values[2*w+1];
      return weights;
    },ct_worker::log_weight,[&](size_t i,char const* name,double value){statistics(i).record(name,value,1.);});
  }
  void measure() {}
  uint64_t completed_sweeps() const {return exchange_ ? exchange_->completed_sweeps() : visit(0,[](auto const& value){return value.completed_sweeps();});}
  double fraction_completed() const {return exchange_ ? exchange_->fraction_completed() : visit(0,[](auto const& value){return value.fraction_completed();});}
  size_t num_sites() const {return sites_;}
  void synchronize() {
    group_.synchronize_physical([&](auto& ar) {
      for (size_t w=0;w<replicas();++w) if (group_.owns_walker(w))
        visit(w,[&](auto const& value){ar["/replicas/"+std::to_string(w)]<<value;});
    },checkpoint_,"/replicas");
  }
  void save(alps::hdf5::archive& ar) const;
  void load(alps::hdf5::archive& ar);

  struct replica_view {
    application const& owner; size_t index;
    size_t num_sites() const {return owner.num_sites();}
    double density_reference() const {return std::numeric_limits<double>::quiet_NaN();}
    auto sampling_parameters() const {
      auto p=owner.parameters_;
      if (owner.exchange_) {p["T"]=1/owner.exchange_->beta(index);p.erase("BETA");}
      return p;
    }
    uint64_t completed_sweeps() const {return owner.completed_sweeps();}
    double fraction_completed() const {return owner.fraction_completed();}
    auto const& get_measurements() const {return owner.statistics(index).get_measurements();}
    auto const& signed_measurements() const {return owner.statistics(index).signed_measurements();}
    auto const& measurement_labels() const {return owner.statistics(index).measurement_labels();}
    template<class T=alps::mc::batch> auto measurement(std::string const& name) const {return owner.statistics(index).template measurement<T>(name);}
    template<class T> auto collect_results_as(alps::mcbase::result_names_type const& names) const {return owner.statistics(index).template collect_results_as<T>(names);}
  };
};
}

namespace looper {
inline void application::save(alps::hdf5::archive& ar) const {
  if (!exchange_) {visit(0,[&](auto const& value){value.save(ar);});return;}
  exchange_->save(ar);
  auto context=ar.get_context();
  for (size_t i=0;i<replicas();++i) {
    auto path="replicas/"+std::to_string(i);
    if (group_.distributed) checkpoint_->copy("/"+path+"/checkpoint",ar,path+"/checkpoint");
    else visit(i,[&](auto const& value){ar[path] << value;});
    ar.set_context(context+"/"+path);
    statistics(i).save_measurements(ar);
    ar.set_context(context);
  }
  ar["/parameters"] << parameters_;
}
inline void application::load(alps::hdf5::archive& ar) {
  application restored(parameters_,bins_,chain_,group_);
  if (!exchange_) restored.visit(0,[&](auto& value){value.load(ar);});
  else {
    auto& r=restored;
    r.exchange_->load(ar);
    if (ar.list_children("replicas").size()!=replicas())
      throw std::invalid_argument("Invalid replica checkpoint state");
    for (size_t i=0;i<replicas();++i) {
      // Expansion order is nonnegative; classical energies in the shared
      // exchange state can have either sign.
      if (r.exchange_->weight_sums()[i][0]<0)
        throw std::invalid_argument("Invalid replica checkpoint feedback");
      if (group_.owns_walker(i)) r.visit(i,[&](auto& value){
        ar["replicas/"+std::to_string(i)] >> value;
        if (value.completed_sweeps()!=r.exchange_->completed_sweeps() || value.measurement("Temperature")->count()!=r.exchange_->production_sweeps())
          throw std::invalid_argument("Inconsistent replica checkpoint sweep counts");
        value.get_measurements().clear();
      });
      auto context=ar.get_context();ar.set_context(context+"/replicas/"+std::to_string(i));
      r.statistics(i).validate_measurements(ar,r.exchange_->completed_sweeps());
      r.statistics(i).alps::mcbase::load(ar);
      ar.set_context(context);
      if (r.statistics(i).measurement("Temperature")->count()!=r.exchange_->production_sweeps())
        throw std::invalid_argument("Inconsistent replica checkpoint sweep counts");
    }
  }
  *this=std::move(restored);
}
}
