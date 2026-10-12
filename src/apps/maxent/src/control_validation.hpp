// SPDX-License-Identifier: MIT
#pragma once
#include <alps/ngs/params.hpp>
#include <cmath>
#include <limits>
#include <stdexcept>
namespace maxent_detail {
inline const alps::params& validate_controls(const alps::params& p) {
    const int n = p["N_ALPHA"], iterations = p["MAX_IT"] | 1000;
    const double lo = p["ALPHA_MIN"], hi = p["ALPHA_MAX"];
    if (n < 2 || iterations < 1 || !std::isfinite(lo) || !std::isfinite(hi) || !(hi > lo && lo > 0))
        throw std::invalid_argument("require N_ALPHA >= 2, MAX_IT >= 1 and 0 < ALPHA_MIN < ALPHA_MAX");
    if (p.defined("T") && p.defined("BETA")) {
        const double t = p["T"], beta = p["BETA"];
        if (!std::isfinite(t*beta) || std::abs(t*beta - 1) > 64*std::numeric_limits<double>::epsilon())
            throw std::invalid_argument("T and BETA must be reciprocal");
    }
    if (p["ENFORCE_NORMALIZATION"] | false) {
        const double sigma = p["SIGMA_NORMALIZATION"];
        if (!std::isfinite(sigma) || sigma <= 0)
            throw std::invalid_argument("SIGMA_NORMALIZATION must be finite and positive");
    }
    return p;
}
}
