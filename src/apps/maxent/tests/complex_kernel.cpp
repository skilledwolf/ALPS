// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#include "maxent_parms.hpp"
#include <boost/math/constants/constants.hpp>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace {
struct kernel_probe : ContiParameters {
    using ContiParameters::ContiParameters;
    void initialize(const alps::params& parameters, const vector_type& frequencies) {
        // Make an unwritten matrix element fail deterministically, regardless of
        // what a previous solver invocation left in the allocator's storage.
        K_ = boost::numeric::ublas::scalar_matrix<double>(
            ndat(), frequencies.size(), std::numeric_limits<double>::quiet_NaN());
        setup_kernel(parameters, frequencies.size(), frequencies);
    }
};
}

int main() {
    const double beta = 2.;
    const double pi = boost::math::constants::pi<double>();
    ContiParameters::vector_type frequencies(4);
    frequencies[0] = -2.;
    frequencies[1] = -0.5;
    frequencies[2] = 0.5;
    frequencies[3] = 2.;
    const alps::maxent::data data{{0., 0., 0., 0.}, {1., 1., 1., 1.}};
    for (const std::string kind : {"fermionic", "bosonic", "anomalous"}) {
        alps::params supplied;
        supplied["BETA"] = beta;
        supplied["NFREQ"] = 4;
        supplied["OMEGA_MAX"] = 3.;
        supplied["DATASPACE"] = "frequency";
        supplied["KERNEL"] = kind;
        supplied["PARTICLE_HOLE_SYMMETRY"] = false;
        const auto parameters = alps::maxent::prepare(supplied, data);
        kernel_probe probe(parameters, data);
        probe.initialize(parameters, frequencies);
        for (int i = 0; i < 2; ++i) {
            const std::complex<double> iw(0., (2 * i + (kind != "bosonic")) * pi / beta);
            for (int j = 0; j < 4; ++j) {
                const double omega = frequencies[j];
                // Current ALPS conventions: G=1/(iw-w), chi=w/(iw-w),
                // anomalous=-w/(iw-w). Review any change during replacement.
                const auto expected = (kind == "fermionic" ? 1. :
                                       kind == "bosonic" ? omega : -omega) / (iw - omega);
                const std::complex<double> actual(probe.K(2 * i, j), probe.K(2 * i + 1, j));
                if (!std::isfinite(actual.real()) || !std::isfinite(actual.imag()) ||
                    std::abs(actual - expected) > 1e-13)
                    throw std::runtime_error(kind + " kernel disagrees at frequency column " +
                                             std::to_string(j));
            }
        }
    }
}
