// Copyright (C) 1997-2010 Synge Todo; 2026 ALPS Collaboration.
// SPDX-License-Identifier: MIT
#pragma once
#include <alps/params.hpp>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace alps::mc {
// Scientific temperature feedback formerly embedded in Parapack's scheduler.
class temperature_grid {
  std::vector<double> beta_;
  static void validate(std::vector<double> const& beta,bool ordered=true) {
    if (beta.size()<2) throw std::invalid_argument("Replica exchange requires at least two temperatures");
    for (size_t i=0;i<beta.size();++i)
      if (!std::isfinite(beta[i]) || beta[i]<=0 || (ordered && i && beta[i]<=beta[i-1]))
        throw std::invalid_argument("Replica inverse temperatures must be positive, finite and distinct");
  }
public:
  explicit temperature_grid(std::vector<double> beta):beta_(std::move(beta)) {
    validate(beta_,false); std::sort(beta_.begin(),beta_.end()); validate(beta_);
  }
  explicit temperature_grid(alps::params const& p) {
    if (p.exists("INVERSE_TEMPERATURE_SET")) beta_=p["INVERSE_TEMPERATURE_SET"].as<std::vector<double>>();
    else if (p.exists("TEMPERATURE_SET")) {
      beta_=p["TEMPERATURE_SET"].as<std::vector<double>>();
      for (auto& value:beta_) value=1/value;
    } else {
      auto n=p["NUM_REPLICAS"].as<size_t>();
      if (n<2) throw std::invalid_argument("NUM_REPLICAS must be at least two");
      bool inverse=p.exists("BETA_MIN") && p.exists("BETA_MAX");
      double low=inverse ? p["BETA_MIN"].as<double>() : 1/p["T_MAX"].as<double>();
      double high=inverse ? p["BETA_MAX"].as<double>() : 1/p["T_MIN"].as<double>();
      if (!std::isfinite(low) || !std::isfinite(high) || low<=0 || high<=low)
        throw std::invalid_argument("Invalid replica temperature endpoints");
      int distribution=p.value_or<int>("TEMPERATURE_DISTRIBUTION_TYPE",inverse ? 1 : 2);
      if (distribution<1 || distribution>3) throw std::invalid_argument("Unknown temperature distribution");
      auto map=[&](double x) { return distribution==2 ? 1/x : distribution==3 ? std::sqrt(x) : x; };
      low=map(low); high=map(high);
      beta_.resize(n);
      for (size_t i=0;i<n;++i) {
        double x=low+(high-low)*double(i)/(n-1);
        beta_[i]=distribution==2 ? 1/x : distribution==3 ? x*x : x;
      }
    }
    validate(beta_,false); std::sort(beta_.begin(),beta_.end()); validate(beta_);
  }
  auto const& values() const { return beta_; }
  size_t size() const { return beta_.size(); }
  double operator[](size_t i) const { return beta_.at(i); }
  void restore(std::vector<double> beta) {
    validate(beta);
    if (beta.size()!=size() || beta.front()!=beta_.front() || beta.back()!=beta_.back())
      throw std::invalid_argument("Checkpoint replica temperature endpoints changed");
    beta_=std::move(beta);
  }

  // Hukushima's neighboring-rate feedback, PRE 60, 3606 (1999). The weight
  // may include field energy as well as expansion order; it must be linear
  // under interpolation, while log_weight supplies the beta dependence.
  template<class Weight,class LogWeight>
  void optimize_rate(std::vector<Weight> const& weights,LogWeight log_weight) {
    if (weights.size()!=size()) throw std::invalid_argument("Replica weight count mismatch");
    auto old=beta_;
    auto interpolate=[&](double beta)->Weight {
      size_t i=std::lower_bound(old.begin(),old.end(),beta)-old.begin();
      i=std::clamp(i,size_t(1),size()-1);
      double t=(beta-old[i-1])/(old[i]-old[i-1]);
      return (1-t)*weights[i-1]+t*weights[i];
    };
    for (int pass=0;pass<64;++pass) for (size_t i=1+pass%2;i+1<size();i+=2) {
      Weight left=interpolate(beta_[i-1]),right=interpolate(beta_[i+1]);
      double low=beta_[i-1],middle=beta_[i],high=beta_[i+1],tolerance=.01*(high-low);
      Weight center=interpolate(middle);
      if (log_weight(left,beta_[i-1])>=log_weight(center,beta_[i-1]) ||
          log_weight(center,beta_[i+1])>=log_weight(right,beta_[i+1])) continue;
      for (int iteration=0;iteration<64;++iteration) {
        double a=log_weight(left,beta_[i-1])+log_weight(center,middle)-log_weight(left,middle)-log_weight(center,beta_[i-1]);
        double b=log_weight(center,middle)+log_weight(right,beta_[i+1])-log_weight(center,beta_[i+1])-log_weight(right,middle);
        if (a<b) { low=middle; middle=(middle+high)/2; }
        else { high=middle; middle=(low+middle)/2; }
        if (high-low<tolerance) break;
        center=interpolate(middle);
      }
      beta_[i]=middle;
    }
    validate(beta_);
  }

  // Invert the integrated feedback density in temperature coordinates. This
  // is the original three-point-regression population method, expressed as a
  // piecewise constant density instead of mutating interval remainders.
  bool optimize_population(std::vector<double> const& upward) {
    if (size()<3 || upward.size()!=size()) throw std::invalid_argument("Population optimization requires at least three replicas");
    for (size_t i=0;i<size();++i)
      if (!std::isfinite(upward[i]) || upward[i]<0 || upward[i]>1 || (i && upward[i]<=upward[i-1])) return false;
    size_t n=size()-1;
    std::vector<double> temperature(size()),fraction(size()),mass(n),optimized(size());
    for (size_t i=0;i<size();++i) { temperature[i]=1/beta_[n-i]; fraction[i]=upward[n-i]; }
    double total=0;
    for (size_t i=0;i<n;++i) {
      size_t first=i ? i-1 : 0;
      double x=0,y=0,xx=0,xy=0;
      for (size_t j=first;j<first+3;++j) { x+=temperature[j]/3; y+=fraction[j]/3; }
      for (size_t j=first;j<first+3;++j) {
        double dx=temperature[j]-x; xx+=dx*dx; xy+=dx*(fraction[j]-y);
      }
      double slope=xy/xx;
      if (!std::isfinite(slope) || slope>=0) return false;
      mass[i]=std::sqrt(-slope*(temperature[i+1]-temperature[i]));
      total+=mass[i];
    }
    if (!std::isfinite(total) || total<=0) return false;
    optimized.front()=beta_.front(); optimized.back()=beta_.back();
    size_t interval=0;
    double consumed=0;
    for (size_t i=1;i<n;++i) {
      double target=total*double(i)/n;
      while (interval+1<n && consumed+mass[interval]<target) consumed+=mass[interval++];
      double t=temperature[interval]+(temperature[interval+1]-temperature[interval])*(target-consumed)/mass[interval];
      optimized[n-i]=1/t;
    }
    validate(optimized); beta_=std::move(optimized);
    return true;
  }
};
}
