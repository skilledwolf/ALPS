// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#pragma once
#include <alps/solvers.hpp>
#include <alps/run_config.hpp>
#include <string_view>
#include <vector>
namespace alps::maxent {
// Owning scientific input; covariance is a flattened row-major square matrix.
struct data {
    std::vector<double> values, errors, covariance, tau, prior_omega, prior_density;
};
std::string_view schema();
data read_data(const params &input);
params prepare(const params &supplied, const data &input);
data prepare_run(run_configuration &run);
} // namespace alps::maxent
