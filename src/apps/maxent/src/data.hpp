// SPDX-License-Identifier: MIT
// Private adaptation of the input record in ALPS reconciliation 098cf7f61
// (src/alps/solvers/include/alps/maxent.hpp). Keep the public params API.
#pragma once
#include <alps/ngs/params.hpp>
#include <vector>
namespace maxent_detail {
struct data {
    std::vector<double> values, errors;
    // Square row-major covariance in input units; used instead of errors.
    std::vector<double> covariance;
    std::vector<double> tau;
    std::vector<double> prior_omega, prior_density;
};
void validate_data(const data&);
data legacy_data(const alps::params&);
int data_count(const data&);
}
