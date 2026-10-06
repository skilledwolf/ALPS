// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#pragma once
#include "loop_worker.h"
#include <alps/mc/replica_exchange.hpp>
#include <optional>
#include <variant>

namespace looper {
class application {
  using ct_worker=loop_worker<path_integral>;
  using sse_worker=loop_worker<sse>;
  using worker=std::variant<std::unique_ptr<ct_worker>,std::unique_ptr<sse_worker>>;
  using weight=ct_worker::weight_parameter_type;
  // Physical walker i and temperature-slot i share storage. Exchange routes
  // each configuration's samples to the statistics at its current temperature.
  std::vector<worker> workers_;
  alps::params parameters_;
  size_t bins_,chain_;
  std::optional<alps::mc::replica_exchange<weight>> exchange_;

  template<class F> decltype(auto) visit(size_t i,F&& f) const {
    return std::visit([&](auto const& value)->decltype(auto){return f(*value);},workers_.at(i));
  }
  native_qmc::simulation& statistics(size_t i) const {
    return visit(i,[](auto& value)->native_qmc::simulation& {return value;});
  }
public:
  application(alps::params const& p,size_t bins,size_t chain):parameters_(p),bins_(bins),chain_(chain) {
    auto algorithm=p.value_or<std::string>("ALGORITHM","loop");
    bool ladder=algorithm.find("exchange")!=std::string::npos;
    if (ladder) exchange_.emplace(p,chain,weight::Zero());
    size_t n=ladder ? exchange_->size() : 1;
    for (size_t i=0;i<n;++i) {
      // Sampling beta is separate from the common XML model parameters.
      double initial_beta=ladder ? exchange_->beta(0) : 0.;
      size_t seed_offset=ladder ? chain*(n+1)+i : chain;
      if (algorithm.find("sse")!=std::string::npos) workers_.push_back(std::make_unique<sse_worker>(p,bins,seed_offset,initial_beta));
      else workers_.push_back(std::make_unique<ct_worker>(p,bins,seed_offset,initial_beta));
    }
    if (ladder) exchange_->init_diagnostics([&](size_t i,char const* name){statistics(i).add_measurement(name,1,false);});
  }
  static alps::params checkpoint_parameters(alps::params p) {return native_qmc::simulation::checkpoint_parameters(p);}
  size_t replicas() const {return workers_.size();}
  bool ladder() const {return bool(exchange_);}
  auto temperatures() const {return exchange_->values();}
  void update() {
    if (!exchange_) {visit(0,[](auto& value){value.update();});return;}
    exchange_->step([&](auto const& walkers,auto const& betas,bool sampling) {
      for (size_t i=0;i<replicas();++i) {
        auto& stats=statistics(i);stats.record_measurements(sampling);
        visit(walkers[i],[&](auto& value){value.run(stats,betas[i]);});
      }
    },[&] {
      std::vector<weight> weights(replicas());
      for (size_t w=0;w<replicas();++w) weights[w]=visit(w,[](auto const& value){return value.weight_parameter();});
      return weights;
    },ct_worker::log_weight,[&](size_t i,char const* name,double value){statistics(i).record(name,value,1.);});
  }
  void measure() {}
  uint64_t completed_sweeps() const {return visit(0,[](auto const& value){return value.completed_sweeps();});}
  double fraction_completed() const {return exchange_ ? exchange_->fraction_completed() : visit(0,[](auto const& value){return value.fraction_completed();});}
  size_t num_sites() const {return visit(0,[](auto const& value){return value.site_count();});}
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
  for (size_t i=0;i<replicas();++i)
    visit(i,[&](auto const& value){ar["replicas/"+std::to_string(i)] << value;});
  ar["/parameters"] << parameters_;
}
inline void application::load(alps::hdf5::archive& ar) {
  application restored(parameters_,bins_,chain_);
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
      r.visit(i,[&](auto& value){
        ar["replicas/"+std::to_string(i)] >> value;
        if (value.completed_sweeps()!=r.exchange_->completed_sweeps() || value.measurement("Temperature")->count()!=r.exchange_->production_sweeps())
          throw std::invalid_argument("Inconsistent replica checkpoint sweep counts");
      });
    }
  }
  *this=std::move(restored);
}
}
