// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#pragma once
#include <alps/mc/physical_moments.hpp>

namespace looper {
using results = std::map<std::string,alps::alea::batch_result<double>>;
// Join aligned raw bins before applying either sign reweighting or nonlinear
// thermodynamics. Joining already jackknifed ratios would propagate bias twice.
inline void derive(results const& raw,results& output,alps::mc::unavailable_results& unavailable) {
  auto ratio=[](auto const& x,size_t offset,size_t size) { return x[offset]/(size==2 ? x[offset+1] : 1.); };
  if (raw.count("Energy") && raw.count("Energy^2")) {
    auto const& energy=raw.at("Energy");
    auto const& square=raw.at("Energy^2");
    double beta=raw.at("Inverse Temperature").mean()[0],volume=raw.at("Volume").mean()[0];
    size_t n=energy.size(),m=square.size();
    alps::mc::estimate(output,&unavailable,"Specific Heat",alps::alea::join(energy,square),
      [=](auto const& x) { return beta*beta/volume*(ratio(x,n,m)-std::pow(ratio(x,0,n),2)); });
  }
  for (std::string name:{"Magnetization","Staggered Magnetization","Generalized Magnetization","Generalized Staggered Magnetization"}) {
    if (!raw.count(name+"^2") || !raw.count(name+"^4")) continue;
    auto const& second=raw.at(name+"^2");
    auto const& fourth=raw.at(name+"^4");
    size_t n=second.size(),m=fourth.size();
    alps::mc::estimate(output,&unavailable,"Binder Ratio of "+name,alps::alea::join(second,fourth),
      [=](auto const& x) { return std::pow(ratio(x,0,n),2)/ratio(x,n,m); });
  }
}
}
