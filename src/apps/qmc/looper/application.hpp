// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#pragma once
#include "loop_worker.h"
#include "temperature_grid.hpp"
#include <numeric>
#include <optional>
#include <variant>

namespace looper {
class application {
  using ct_worker=loop_worker<path_integral>;
  using sse_worker=loop_worker<sse>;
  using worker=std::variant<std::unique_ptr<ct_worker>,std::unique_ptr<sse_worker>>;
  using weight=ct_worker::weight_parameter_type;
  // Physical walker i and temperature-slot i share storage, but updates use
  // walker_at_ to route each configuration's samples to its current temperature.
  std::vector<worker> workers_;
  alps::params parameters_;
  size_t bins_,chain_;
  std::optional<temperature_grid> grid_;
  alps::random01 random_;
  bool exchange_=false,random_exchange_=false,population_=false,ready_=false;
  uint64_t steps_=0,production_=0,warm_=0,interval_=1,stage_=0,stage_count_=0,events_=0,returnees_=0;
  double factor_=1;
  std::vector<uint64_t> blocks_;
  std::vector<size_t> walker_at_; // temperature -> physical walker
  std::vector<int> direction_; // 0: unlabelled, 1: downward, 2: upward
  std::vector<weight> weight_sum_;
  std::vector<double> up_,down_;

  template<class F> decltype(auto) visit(size_t i,F&& f) const {
    return std::visit([&](auto const& value)->decltype(auto){return f(*value);},workers_.at(i));
  }
  bool optimizing() const { return stage_<blocks_.size(); }
  native_qmc::simulation& statistics(size_t i) const {
    return visit(i,[](auto& value)->native_qmc::simulation& {return value;});
  }
  uint64_t grow(uint64_t count) const {
    long double next=std::ceil(std::max(1.,factor_)*count);
    if (!std::isfinite(next) || next>=std::ldexp(1.L,64) || count==std::numeric_limits<uint64_t>::max())
      throw std::overflow_error("Replica optimization sweep count overflow");
    return std::max(count+1,uint64_t(next));
  }
  // Stage zero equilibrates the initial ladder; subsequent stages move its
  // interior temperatures. Production begins only after feedback and warmup.
  void feedback() {
    if (!optimizing() || stage_count_<blocks_[stage_]) return;
    bool success=true;
    if (population_) {
      success=returnees_>=replicas() && std::find(direction_.begin(),direction_.end(),0)==direction_.end();
      std::vector<double> fraction(replicas());
      for (size_t i=0;i<replicas();++i) fraction[i]=up_[i]/(up_[i]+down_[i]);
      if (stage_ && success) success=grid_->optimize_population(fraction);
    } else if (stage_) {
      auto means=weight_sum_;
      for (auto& value:means) value/=double(events_);
      grid_->optimize_rate(means,ct_worker::log_weight);
    }
    if (!success) { for (size_t i=stage_;i<blocks_.size();++i) blocks_[i]=grow(blocks_[i]); return; }
    ++stage_; stage_count_=0; events_=0; returnees_=0;
    std::fill(weight_sum_.begin(),weight_sum_.end(),weight::Zero());
    std::fill(up_.begin(),up_.end(),0.); std::fill(down_.begin(),down_.end(),0.);
    if (!optimizing() && !parameters_["THERMALIZATION"].as<uint64_t>()) ready_=true;
  }
  void exchange() {
    std::vector<weight> weights(replicas());
    for (size_t w=0;w<replicas();++w) weights[w]=visit(w,[](auto const& value){return value.weight_parameter();});
    if (optimizing()) {
      ++events_;
      for (size_t i=0;i<replicas();++i) weight_sum_[i]+=weights[walker_at_[i]];
    }
    std::vector<size_t> edges;
    if (random_exchange_) {
      edges.resize(replicas()-1); std::iota(edges.begin(),edges.end(),0);
      for (size_t i=edges.size();i>1;--i) std::swap(edges[i-1],edges[size_t(random_()*i)]);
    } else for (size_t i=(steps_/interval_)%2;i+1<replicas();i+=2) edges.push_back(i);
    for (size_t i:edges) {
      auto a=walker_at_[i],b=walker_at_[i+1];
      double logp=ct_worker::log_weight(weights[b]-weights[a],(*grid_)[i])-
                  ct_worker::log_weight(weights[b]-weights[a],(*grid_)[i+1]);
      bool accepted=std::log(random_())<logp;
      if (accepted) std::swap(walker_at_[i],walker_at_[i+1]);
      statistics(i).record("EXMC: Acceptance Rate",double(accepted),1.);
    }
    auto first=walker_at_.front(),last=walker_at_.back();
    bool returned=direction_[first]==2;
    for (size_t w=0;w<replicas();++w)
      statistics(w).record("EXMC: Inverse Round-Trip Time",double(returned && w==first),1.);
    statistics(0).record("EXMC: Average Inverse Round-Trip Time",double(returned)/replicas(),1.);
    returnees_+=returned;
    direction_[first]=1;
    if (direction_[last]==1) direction_[last]=2;
    for (size_t i=0;i<replicas();++i) {
      double up=direction_[walker_at_[i]]==2,down=direction_[walker_at_[i]]==1;
      up_[i]+=up; down_[i]+=down;
      statistics(i).record("EXMC: Ratio of Upward-Moving Walker",up,1.);
      statistics(i).record("EXMC: Ratio of Downward-Moving Walker",down,1.);
    }
    feedback();
  }
public:
  application(alps::params const& p,size_t bins,size_t chain):parameters_(p),bins_(bins),chain_(chain),random_(0,p.value_or<std::string>("RNG","mt19937")) {
    auto algorithm=p.value_or<std::string>("ALGORITHM","loop");
    bool ladder=algorithm.find("exchange")!=std::string::npos;
    if (ladder) grid_.emplace(p);
    size_t n=grid_ ? grid_->size() : 1;
    if (ladder) {
      auto seed=p.value_or<uint64_t>("SEED",42);
      if (n>INT_MAX-seed || chain>(INT_MAX-seed-n)/(n+1))
        throw std::invalid_argument("Replica RNG seed range exceeded");
      random_.seed(uint32_t(seed+chain*(n+1)+n));
      exchange_=!p.value_or("NO_EXCHANGE",false);
      random_exchange_=p.value_or("RANDOM_EXCHANGE",false);
      interval_=p.value_or<uint64_t>("EXCHANGE_INTERVAL",1);
      if (!interval_) throw std::invalid_argument("EXCHANGE_INTERVAL must be positive");
      bool optimize=p.value_or("OPTIMIZE_TEMPERATURE",false) || p.value_or("TEMPERATURE_OPTIMIZATION",false);
      if (optimize) {
        if (!exchange_) throw std::invalid_argument("Temperature optimization requires replica exchange");
        auto type=p.value_or<std::string>("OPTIMIZATION_TYPE","rate");
        population_=type=="population";
        if (type!="rate" && !population_) throw std::invalid_argument("Unknown temperature optimization type");
        if (population_ && n<3) throw std::invalid_argument("Population optimization requires at least three replicas");
        factor_=p.value_or<double>("BLOCK_SWEEP_FACTOR",population_ ? 2. : 1.);
        if (!std::isfinite(factor_) || factor_<1) throw std::invalid_argument("BLOCK_SWEEP_FACTOR must be finite and at least one");
        auto iterations=p.value_or<size_t>("OPTIMIZATION_ITERATIONS",population_ ? 7 : 1);
        uint64_t block=p.value_or<uint64_t>("INITIAL_BLOCK_SWEEPS",population_ ? 512 : std::max(uint64_t(1),p["THERMALIZATION"].as<uint64_t>()));
        if (!block) throw std::invalid_argument("INITIAL_BLOCK_SWEEPS must be positive");
        for (size_t i=0;i<=iterations;++i) {
          blocks_.push_back(block);
          if (i<iterations && factor_>1) block=grow(block);
        }
      }
    }
    ready_=!optimizing() && !p["THERMALIZATION"].as<uint64_t>();
    walker_at_.resize(n); std::iota(walker_at_.begin(),walker_at_.end(),0);
    direction_.assign(n,0); direction_[0]=1;
    weight_sum_.assign(n,weight::Zero()); up_.assign(n,0.); down_.assign(n,0.);
    for (size_t i=0;i<n;++i) {
      // Sampling beta is separate from the common XML model parameters.
      double initial_beta=ladder ? (*grid_)[0] : 0.;
      size_t seed_offset=ladder ? chain*(n+1)+i : chain;
      if (algorithm.find("sse")!=std::string::npos) workers_.push_back(std::make_unique<sse_worker>(p,bins,seed_offset,initial_beta));
      else workers_.push_back(std::make_unique<ct_worker>(p,bins,seed_offset,initial_beta));
      if (ladder) {
        auto& stats=statistics(i);
        for (auto name:{"EXMC: Temperature","EXMC: Inverse Temperature"}) stats.add_measurement(name,1,false);
        if (exchange_) {
          for (auto name:{"EXMC: Ratio of Upward-Moving Walker","EXMC: Ratio of Downward-Moving Walker","EXMC: Inverse Round-Trip Time"}) stats.add_measurement(name,1,false);
          if (i+1<n) stats.add_measurement("EXMC: Acceptance Rate",1,false);
          if (!i) stats.add_measurement("EXMC: Average Inverse Round-Trip Time",1,false);
        }
      }
    }
  }
  static alps::params checkpoint_parameters(alps::params p) {return native_qmc::simulation::checkpoint_parameters(p);}
  size_t replicas() const {return workers_.size();}
  bool ladder() const {return bool(grid_);}
  auto temperatures() const {return grid_->values();}
  void update() {
    if (!grid_) {visit(0,[](auto& value){value.update();});return;}
    bool sampling=ready_,was_optimizing=optimizing();
    ++steps_; if (was_optimizing) ++stage_count_;
    for (size_t i=0;i<replicas();++i) {
      auto& stats=statistics(i);stats.record_measurements(sampling);
      visit(walker_at_[i],[&](auto& value){value.run(stats,(*grid_)[i]);});
      stats.record("EXMC: Temperature",1/(*grid_)[i],1.);
      stats.record("EXMC: Inverse Temperature",(*grid_)[i],1.);
    }
    if (exchange_ && steps_%interval_==0) exchange();
    if (sampling) ++production_;
    else if (!was_optimizing) {
      ++warm_;
      if (warm_>=parameters_["THERMALIZATION"].as<uint64_t>()) ready_=true;
    }
  }
  void measure() {}
  uint64_t completed_sweeps() const {return visit(0,[](auto const& value){return value.completed_sweeps();});}
  double fraction_completed() const {return grid_ ? double(production_)/parameters_["SWEEPS"].as<uint64_t>() : visit(0,[](auto const& value){return value.fraction_completed();});}
  size_t num_sites() const {return visit(0,[](auto const& value){return value.site_count();});}
  void save(alps::hdf5::archive& ar) const;
  void load(alps::hdf5::archive& ar);

  struct replica_view {
    application const& owner; size_t index;
    size_t num_sites() const {return owner.num_sites();}
    double density_reference() const {return std::numeric_limits<double>::quiet_NaN();}
    auto sampling_parameters() const {
      auto p=owner.parameters_;
      if (owner.grid_) {p["T"]=1/(*owner.grid_)[index];p.erase("BETA");}
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
  if (!grid_) {visit(0,[&](auto const& value){value.save(ar);});return;}
  ar["exchange/version"] << uint64_t(1);
  ar["exchange/beta"] << grid_->values();
  ar["exchange/random"] << random_;
  ar["exchange/steps"] << steps_; ar["exchange/production"] << production_;
  ar["exchange/warmup"] << warm_; ar["exchange/ready"] << ready_;
  ar["exchange/stage"] << stage_; ar["exchange/stage_count"] << stage_count_;
  ar["exchange/blocks"] << blocks_; ar["exchange/events"] << events_;
  ar["exchange/returnees"] << returnees_; ar["exchange/walkers"] << walker_at_;
  ar["exchange/direction"] << direction_; ar["exchange/up"] << up_; ar["exchange/down"] << down_;
  std::vector<std::array<double,2>> weights;
  for (auto const& value:weight_sum_) weights.push_back({value[0],value[1]});
  ar["exchange/weights"] << weights;
  for (size_t i=0;i<replicas();++i)
    visit(i,[&](auto const& value){ar["replicas/"+std::to_string(i)] << value;});
  ar["/parameters"] << parameters_;
}
inline void application::load(alps::hdf5::archive& ar) {
  application restored(parameters_,bins_,chain_);
  if (!grid_) restored.visit(0,[&](auto& value){value.load(ar);});
  else {
    auto& r=restored;
    uint64_t version;
    std::vector<double> beta;
    std::vector<std::array<double,2>> weights;
    ar["exchange/version"] >> version;
    ar["exchange/beta"] >> beta; r.grid_->restore(std::move(beta));
    ar["exchange/random"] >> r.random_;
    ar["exchange/steps"] >> r.steps_; ar["exchange/production"] >> r.production_;
    ar["exchange/warmup"] >> r.warm_; ar["exchange/ready"] >> r.ready_;
    ar["exchange/stage"] >> r.stage_; ar["exchange/stage_count"] >> r.stage_count_;
    ar["exchange/blocks"] >> r.blocks_; ar["exchange/events"] >> r.events_;
    ar["exchange/returnees"] >> r.returnees_; ar["exchange/walkers"] >> r.walker_at_;
    ar["exchange/direction"] >> r.direction_; ar["exchange/up"] >> r.up_; ar["exchange/down"] >> r.down_;
    ar["exchange/weights"] >> weights;
    if (version!=1 || r.random_.name()!=random_.name() || r.stage_>r.blocks_.size() ||
        r.blocks_.size()!=blocks_.size() || r.production_>parameters_["SWEEPS"].as<uint64_t>() ||
        r.warm_>r.steps_ || r.production_>r.steps_-r.warm_ || r.stage_count_>r.steps_ ||
        r.events_>r.steps_/interval_ || r.returnees_>r.steps_/interval_ ||
        (!r.ready_ && r.production_) || (r.ready_ && (r.optimizing() || r.warm_<parameters_["THERMALIZATION"].as<uint64_t>())) ||
        r.walker_at_.size()!=replicas() || r.direction_.size()!=replicas() ||
        r.up_.size()!=replicas() || r.down_.size()!=replicas() || weights.size()!=replicas() ||
        ar.list_children("replicas").size()!=replicas())
      throw std::invalid_argument("Invalid replica checkpoint state");
    for (size_t i=0;i<blocks_.size();++i)
      if (r.blocks_[i]<blocks_[i]) throw std::invalid_argument("Invalid replica optimization schedule");
    auto permutation=r.walker_at_; std::sort(permutation.begin(),permutation.end());
    for (size_t i=0;i<replicas();++i) {
      if (permutation[i]!=i || r.direction_[i]<0 || r.direction_[i]>2 ||
          !std::isfinite(r.up_[i]) || !std::isfinite(r.down_[i]) || r.up_[i]<0 || r.down_[i]<0 ||
          r.up_[i]+r.down_[i]>double(r.steps_/interval_) ||
          !std::isfinite(weights[i][0]) || !std::isfinite(weights[i][1]) || weights[i][0]<0)
        throw std::invalid_argument("Invalid replica checkpoint feedback");
      r.weight_sum_[i]={weights[i][0],weights[i][1]};
      r.visit(i,[&](auto& value){
        ar["replicas/"+std::to_string(i)] >> value;
        if (value.completed_sweeps()!=r.steps_ || value.measurement("Temperature")->count()!=r.production_)
          throw std::invalid_argument("Inconsistent replica checkpoint sweep counts");
      });
    }
  }
  *this=std::move(restored);
}
}
