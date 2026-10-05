// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#pragma once
#include "../mc/physical_moments.hpp"
#include <alps/params.hpp>
#include <alps/ngs/make_deprecated_parameters.hpp>
#include <alps/expression.h>
namespace native_qmc {
// Each signed measurement carries its own aligned denominator, including when
// particle-number restrictions or measurement intervals select different samples.
struct divide_sign : alps::alea::transformer<double> {
    size_t size;
    explicit divide_sign(size_t n):size(n){}
    size_t in_size() const override { return size; }
    size_t out_size() const override { return size-1; }
    alps::alea::column<double> operator()(alps::alea::column<double> const& x) const override {
        return x.head(size-1)/x[size-1];
    }
};

inline double inverse_temperature(alps::params const& p) {
    auto legacy=alps::make_deprecated_parameters(p);
    for (auto key:{"Beta","beta","BETA"}) if (legacy.defined(key)) return alps::evaluate<double>(legacy[key],legacy);
    for (auto key:{"T","TEMPERATURE","temperature"}) if (legacy.defined(key)) return 1/alps::evaluate<double>(legacy[key],legacy);
    throw std::invalid_argument("Temperature is missing");
}
inline void derive(std::map<std::string,alps::alea::batch_result<double>>& results,
                   alps::params const& p,size_t sites,native_mc::unavailable_results& unavailable) {
    double beta=inverse_temperature(p);
    if (results.count("Centered Density Moments"))
        native_mc::estimate(results,&unavailable,"Compressibility",results.at("Centered Density Moments"),
            [=](auto const& x){return beta*sites*(x[1]/x[2]-std::pow(x[0]/x[2],2));});
    if (p.value_or("USE_1D_STIFFNESS",false) && results.count("Winding number histogram"))
        native_mc::estimate(results,&unavailable,"Superfluid stiffness (1D estimator)",results.at("Winding number histogram"),
            [=](auto const& x){return p.value_or<double>("L",sites)/(2*beta*std::log(2*x[1]/(x[0]+x[2])));});
}

}
