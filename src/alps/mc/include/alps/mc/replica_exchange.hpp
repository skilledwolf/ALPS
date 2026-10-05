// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#pragma once
#include "temperature_grid.hpp"
#include <alps/hdf5/vector.hpp>
#include <alps/ngs/random01.hpp>
#include <climits>
#include <limits>
#include <numeric>
#include <type_traits>
#include <utility>

namespace alps::mc {
// Exchange owns the temperature assignment and feedback, while applications
// own physical walkers and temperature-slot statistics. Weight is either a
// floating-point scalar energy or a fixed Eigen vector of floating-point
// linear weight parameters. Pass an explicit zero because an Eigen vector's
// default constructor does not initialize its data.
template<class Weight> class replica_exchange {
  static_assert([] {
    if constexpr (std::is_arithmetic_v<Weight>) return std::is_floating_point_v<Weight>;
    else return Weight::IsVectorAtCompileTime && Weight::SizeAtCompileTime>0 && std::is_floating_point_v<typename Weight::Scalar>;
  }(),"Replica exchange requires a floating-point scalar or fixed floating-point Eigen vector");
  temperature_grid grid_;
  alps::random01 random_;
  Weight zero_;
  bool enabled_,random_exchange_,population_=false,ready_=false;
  uint64_t target_,thermalization_,steps_=0,production_=0,warm_=0,interval_,stage_=0,stage_count_=0,events_=0,returnees_=0,initial_block_=0;
  double factor_=1;
  std::vector<uint64_t> blocks_;
  std::vector<size_t> walker_at_; // temperature slot -> physical walker
  std::vector<int> direction_; // 0: unlabelled, 1: downward, 2: upward
  std::vector<Weight> weight_sum_;
  std::vector<double> up_,down_;

  bool optimizing() const {return stage_<blocks_.size();}
  uint64_t grow(uint64_t count) const {
    long double next=std::ceil(std::max(1.,factor_)*count);
    if (!std::isfinite(next) || next>=std::ldexp(1.L,64) || count==std::numeric_limits<uint64_t>::max())
      throw std::overflow_error("Replica optimization sweep count overflow");
    return std::max(count+1,uint64_t(next));
  }
  template<class LogWeight> void feedback(LogWeight const& log_weight) {
    if (!optimizing() || stage_count_<blocks_[stage_]) return;
    bool success=true;
    if (population_) {
      success=returnees_>=size() && std::find(direction_.begin(),direction_.end(),0)==direction_.end();
      std::vector<double> fraction(size());
      for (size_t i=0;i<size();++i) fraction[i]=up_[i]/(up_[i]+down_[i]);
      if (stage_ && success) success=grid_.optimize_population(fraction);
    } else if (stage_) {
      auto means=weight_sum_;
      for (auto& value:means) value/=double(events_);
      grid_.optimize_rate(means,log_weight);
    }
    if (!success) {for (size_t i=stage_;i<blocks_.size();++i) blocks_[i]=grow(blocks_[i]);return;}
    ++stage_;stage_count_=0;events_=0;returnees_=0;
    std::fill(weight_sum_.begin(),weight_sum_.end(),zero_);
    std::fill(up_.begin(),up_.end(),0.);std::fill(down_.begin(),down_.end(),0.);
    if (!optimizing() && !thermalization_) ready_=true;
  }
  template<class LogWeight,class Record>
  void exchange(std::vector<Weight> const& weights,LogWeight const& log_weight,Record const& record) {
    if (weights.size()!=size()) throw std::invalid_argument("Replica weight count mismatch");
    if (optimizing()) {
      ++events_;
      for (size_t i=0;i<size();++i) weight_sum_[i]+=weights[walker_at_[i]];
    }
    std::vector<size_t> edges;
    if (random_exchange_) {
      edges.resize(size()-1);std::iota(edges.begin(),edges.end(),0);
      for (size_t i=edges.size();i>1;--i) std::swap(edges[i-1],edges[size_t(random_()*i)]);
    } else for (size_t i=(steps_/interval_)%2;i+1<size();i+=2) edges.push_back(i);
    for (size_t i:edges) {
      auto a=walker_at_[i],b=walker_at_[i+1];
      double logp=log_weight(weights[b]-weights[a],grid_[i])-log_weight(weights[b]-weights[a],grid_[i+1]);
      bool accepted=std::log(random_())<logp;
      if (accepted) std::swap(walker_at_[i],walker_at_[i+1]);
      record(i,"EXMC: Acceptance Rate",double(accepted));
    }
    auto first=walker_at_.front(),last=walker_at_.back();
    bool returned=direction_[first]==2;
    // Round-trip histories belong to physical walkers; population and
    // acceptance diagnostics belong to temperature slots.
    for (size_t w=0;w<size();++w)
      record(w,"EXMC: Inverse Round-Trip Time",double(returned && w==first));
    record(0,"EXMC: Average Inverse Round-Trip Time",double(returned)/size());
    returnees_+=returned;direction_[first]=1;
    if (direction_[last]==1) direction_[last]=2;
    for (size_t i=0;i<size();++i) {
      double up=direction_[walker_at_[i]]==2,down=direction_[walker_at_[i]]==1;
      up_[i]+=up;down_[i]+=down;
      record(i,"EXMC: Ratio of Upward-Moving Walker",up);
      record(i,"EXMC: Ratio of Downward-Moving Walker",down);
    }
    feedback(log_weight);
  }
  size_t weight_size() const {
    if constexpr (std::is_arithmetic_v<Weight>) return 1;
    else return zero_.size();
  }
public:
  replica_exchange(alps::params const& p,size_t chain,Weight zero):grid_(p),
      random_(0,p.value_or<std::string>("RNG","mt19937")),zero_(std::move(zero)),
      enabled_(!p.value_or("NO_EXCHANGE",false)),random_exchange_(p.value_or("RANDOM_EXCHANGE",false)),
      target_(p["SWEEPS"].as<uint64_t>()),thermalization_(p["THERMALIZATION"].as<uint64_t>()),
      interval_(p.value_or<uint64_t>("EXCHANGE_INTERVAL",1)) {
    const size_t n=size();
    auto seed=p.value_or<uint64_t>("SEED",42);
    if (seed>INT_MAX || n>INT_MAX-seed || chain>(INT_MAX-seed-n)/(n+1))
      throw std::invalid_argument("Replica RNG seed range exceeded");
    random_.seed(uint32_t(seed+chain*(n+1)+n));
    if (!interval_) throw std::invalid_argument("EXCHANGE_INTERVAL must be positive");
    bool optimize=p.value_or("OPTIMIZE_TEMPERATURE",false) || p.value_or("TEMPERATURE_OPTIMIZATION",false);
    if (optimize) {
      if (!enabled_) throw std::invalid_argument("Temperature optimization requires replica exchange");
      auto type=p.value_or<std::string>("OPTIMIZATION_TYPE","rate");
      population_=type=="population";
      if (type!="rate" && !population_) throw std::invalid_argument("Unknown temperature optimization type");
      if (population_ && n<3) throw std::invalid_argument("Population optimization requires at least three replicas");
      factor_=p.value_or<double>("BLOCK_SWEEP_FACTOR",population_ ? 2. : 1.);
      if (!std::isfinite(factor_) || factor_<1) throw std::invalid_argument("BLOCK_SWEEP_FACTOR must be finite and at least one");
      auto iterations=p.value_or<size_t>("OPTIMIZATION_ITERATIONS",population_ ? 7 : 1);
      uint64_t block=p.value_or<uint64_t>("INITIAL_BLOCK_SWEEPS",population_ ? 512 : std::max(uint64_t(1),thermalization_));
      if (!block) throw std::invalid_argument("INITIAL_BLOCK_SWEEPS must be positive");
      initial_block_=block;
      for (size_t i=0;i<=iterations;++i) {
        blocks_.push_back(block);
        if (i<iterations && factor_>1) block=grow(block);
      }
    }
    ready_=!optimizing() && !thermalization_;
    walker_at_.resize(n);std::iota(walker_at_.begin(),walker_at_.end(),0);
    direction_.assign(n,0);direction_[0]=1;
    weight_sum_.assign(n,zero_);up_.assign(n,0.);down_.assign(n,0.);
  }
  size_t size() const {return grid_.size();}
  double beta(size_t i) const {return grid_[i];}
  auto const& values() const {return grid_.values();}
  bool enabled() const {return enabled_;}
  uint64_t completed_sweeps() const {return steps_;}
  uint64_t production_sweeps() const {return production_;}
  double fraction_completed() const {return double(production_)/target_;}
  auto const& weight_sums() const {return weight_sum_;}
  template<class Add> void init_diagnostics(Add const& add) const {
    for (size_t i=0;i<size();++i) {
      for (auto name:{"EXMC: Temperature","EXMC: Inverse Temperature"}) add(i,name);
      if (enabled_) {
        for (auto name:{"EXMC: Ratio of Upward-Moving Walker","EXMC: Ratio of Downward-Moving Walker","EXMC: Inverse Round-Trip Time"}) add(i,name);
        if (i+1<size()) add(i,"EXMC: Acceptance Rate");
        if (!i) add(i,"EXMC: Average Inverse Round-Trip Time");
      }
    }
  }
  // Update(walker,slot,beta,sampling) routes the physical sample to its slot.
  // Weights() returns physical-walker weights only when an exchange is due;
  // applications may gather them from distributed workers in that callback.
  template<class Update,class Weights,class LogWeight,class Record>
  void step(Update const& update,Weights const& weights,LogWeight const& log_weight,Record const& record) {
    bool sampling=ready_,was_optimizing=optimizing();
    ++steps_;if (was_optimizing) ++stage_count_;
    for (size_t i=0;i<size();++i) {
      update(walker_at_[i],i,grid_[i],sampling);
      record(i,"EXMC: Temperature",1/grid_[i]);record(i,"EXMC: Inverse Temperature",grid_[i]);
    }
    if (enabled_ && steps_%interval_==0) exchange(weights(),log_weight,record);
    if (sampling) ++production_;
    else if (!was_optimizing) {++warm_;if (warm_>=thermalization_) ready_=true;}
  }
  void save(alps::hdf5::archive& ar) const {
    ar["exchange/version"]<<uint64_t(1);ar["exchange/beta"]<<grid_.values();ar["exchange/random"]<<random_;
    ar["exchange/steps"]<<steps_;ar["exchange/production"]<<production_;
    ar["exchange/warmup"]<<warm_;ar["exchange/ready"]<<ready_;
    ar["exchange/stage"]<<stage_;ar["exchange/stage_count"]<<stage_count_;
    ar["exchange/blocks"]<<blocks_;ar["exchange/events"]<<events_;
    ar["exchange/returnees"]<<returnees_;ar["exchange/walkers"]<<walker_at_;
    ar["exchange/direction"]<<direction_;ar["exchange/up"]<<up_;ar["exchange/down"]<<down_;
    // The existing loop checkpoint uses a native N x 2 dataset. Scalar
    // energies use the same rectangular representation with one column.
    std::vector<std::vector<double>> weights(size(),std::vector<double>(weight_size()));
    for (size_t i=0;i<size();++i) {
      if constexpr (std::is_arithmetic_v<Weight>) weights[i][0]=weight_sum_[i];
      else for (size_t j=0;j<weight_size();++j) weights[i][j]=weight_sum_[i][j];
    }
    ar["exchange/weights"]<<weights;
  }
  void load(alps::hdf5::archive& ar) {
    auto r=*this;
    uint64_t version;std::vector<double> beta;std::vector<std::vector<double>> weights;
    ar["exchange/version"]>>version;ar["exchange/beta"]>>beta;r.grid_.restore(std::move(beta));
    ar["exchange/random"]>>r.random_;ar["exchange/steps"]>>r.steps_;ar["exchange/production"]>>r.production_;
    ar["exchange/warmup"]>>r.warm_;ar["exchange/ready"]>>r.ready_;
    ar["exchange/stage"]>>r.stage_;ar["exchange/stage_count"]>>r.stage_count_;
    ar["exchange/blocks"]>>r.blocks_;ar["exchange/events"]>>r.events_;
    ar["exchange/returnees"]>>r.returnees_;ar["exchange/walkers"]>>r.walker_at_;
    ar["exchange/direction"]>>r.direction_;ar["exchange/up"]>>r.up_;ar["exchange/down"]>>r.down_;
    ar["exchange/weights"]>>weights;
    if (version!=1 || r.random_.name()!=random_.name() || r.stage_>r.blocks_.size() ||
        r.blocks_.size()!=blocks_.size() || r.production_>target_ ||
        r.warm_>r.steps_ || r.production_>r.steps_-r.warm_ || r.stage_count_>r.steps_ ||
        r.events_>r.steps_/interval_ || r.returnees_>r.steps_/interval_ ||
        (!r.ready_ && r.production_) || (r.ready_ && (r.optimizing() || r.warm_<thermalization_)) ||
        r.walker_at_.size()!=size() || r.direction_.size()!=size() ||
        r.up_.size()!=size() || r.down_.size()!=size() || weights.size()!=size())
      throw std::invalid_argument("Invalid replica checkpoint state");
    // A checkpoint may rewind feedback that has already grown the live
    // schedule. Validate against the configured minima, not that live state.
    uint64_t minimum=initial_block_;
    for (size_t i=0;i<blocks_.size();++i) {
      if (r.blocks_[i]<minimum) throw std::invalid_argument("Invalid replica optimization schedule");
      if (i+1<blocks_.size() && factor_>1) minimum=grow(minimum);
    }
    auto permutation=r.walker_at_;std::sort(permutation.begin(),permutation.end());
    for (size_t i=0;i<size();++i) {
      if (permutation[i]!=i || r.direction_[i]<0 || r.direction_[i]>2 ||
          !std::isfinite(r.up_[i]) || !std::isfinite(r.down_[i]) || r.up_[i]<0 || r.down_[i]<0 ||
          r.up_[i]+r.down_[i]>double(r.steps_/interval_) || weights[i].size()!=weight_size() ||
          !std::all_of(weights[i].begin(),weights[i].end(),[](double x){return std::isfinite(x);}))
        throw std::invalid_argument("Invalid replica checkpoint feedback");
      if constexpr (std::is_arithmetic_v<Weight>) r.weight_sum_[i]=weights[i][0];
      else for (size_t j=0;j<weight_size();++j) r.weight_sum_[i][j]=weights[i][j];
    }
    *this=std::move(r);
  }
};
}
