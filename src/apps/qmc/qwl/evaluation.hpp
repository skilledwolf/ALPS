// Copyright (C) 2004 Stefan Wessel; 2026 ALPS Collaboration.
// SPDX-License-Identifier: MIT
#pragma once
#include <alps/alea/mean.hpp>
#include <map>
#include <string>
#include <cmath>

namespace qwl {
// Log-sum-exp keeps the truncated partition sum representable at low T.
// Nonzero windows have no absolute normalization: only energy/heat are defined.
inline std::map<std::string,double> evaluate(alps::alea::column<double> const& coefficients,
    double offset,size_t sites,unsigned first,double temperature,
    alps::alea::column<double> const& uniform={},alps::alea::column<double> const& staggered={}) {
  if (!sites || !coefficients.size() || !std::isfinite(temperature) || temperature<=0 || !std::isfinite(offset))
    throw std::invalid_argument("Invalid QWL thermodynamic input");
  alps::alea::column<double> logweights(coefficients.size());
  for (Eigen::Index i=0;i<coefficients.size();++i) {
    if (std::isnan(coefficients[i]) || coefficients[i]==INFINITY)
      throw std::invalid_argument("Invalid QWL coefficient");
    logweights[i]=coefficients[i]-(first+i)*std::log(temperature);
  }
  double maximum=logweights.maxCoeff();
  if (!std::isfinite(maximum)) throw std::invalid_argument("QWL has no sampled coefficients");
  alps::alea::column<double> weights=(logweights.array()-maximum).exp().matrix();
  double z=weights.sum(), n=0, variance=0;
  for (Eigen::Index i=0;i<weights.size();++i) n+=(first+i)*weights[i]/z;
  for (Eigen::Index i=0;i<weights.size();++i) variance+=std::pow(first+i-n,2)*weights[i]/z;
  std::map<std::string,double> result{{"Energy Density",(offset-temperature*n)/sites},
                                    {"Specific Heat per Site",(variance-n)/sites}};
  if (!first) {
    double free=offset/sites-temperature*(std::log(z)+maximum)/sites;
    result["Free Energy Density"]=free;
    result["Entropy Density"]=(result.at("Energy Density")-free)/temperature;
  }
  auto magnetic=[&](alps::alea::column<double> const& values) {
    if (values.size()!=weights.size()) throw std::invalid_argument("QWL magnetic coefficient shape mismatch");
    double total=0;
    for (Eigen::Index i=0;i<values.size();++i) {
      if (!weights[i]) continue;
      if (!std::isfinite(values[i]) || values[i]<0)
        throw std::invalid_argument("QWL magnetic coefficient has no samples; increase production sweeps");
      total+=values[i]*weights[i]/z;
    }
    return total/sites;
  };
  if (uniform.size()) {
    result["Uniform Structure Factor per Site"]=magnetic(uniform);
    result["Uniform Susceptibility per Site"]=magnetic(uniform)/temperature;
  }
  if (staggered.size()) result["Staggered Structure Factor per Site"]=magnetic(staggered);
  return result;
}
}
