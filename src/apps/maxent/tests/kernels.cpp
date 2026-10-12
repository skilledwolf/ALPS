// SPDX-License-Identifier: MIT
#include <gtest/gtest.h>
#include "maxent_parms.hpp"
#include <complex>

namespace {
alps::params parameters(const char* kernel, const char* space, int n = 4) {
    alps::params p;
    p["NFREQ"] = 4; p["NDAT"] = n; p["BETA"] = 2.; p["NORM"] = 1.;
    p["OMEGA_MAX"] = 4.; p["FREQUENCY_GRID"] = "linear";
    p["KERNEL"] = kernel; p["DATASPACE"] = space;
    for (int i = 0; i < n; ++i) {
        p["X_" + std::to_string(i)] = 0.; p["SIGMA_" + std::to_string(i)] = 1.;
    }
    return p;
}
struct kernel_probe : ContiParameters {
    kernel_probe(const alps::params& p, std::initializer_list<double> omega) : ContiParameters(p) {
        vector_type f(omega.size()); std::copy(omega.begin(), omega.end(), f.begin());
        setup_kernel(p, f.size(), f);
    }
};
}
TEST(MaxEntKernel, ComplexColumnsMatchAnalyticPoleResponse) {
    const double omega[] = {-1.7, 0., .8};
    for (const auto* type : {"fermionic", "bosonic", "anomalous"}) {
        auto p = parameters(type, "frequency"); kernel_probe k(p, {-1.7, 0., .8});
        const bool bose = std::string(type) == "bosonic";
        for (int i = 0; i < 2; ++i) for (int j = 0; j < 3; ++j) {
            const double wn = (2*i + (bose ? 0 : 1)) * std::acos(-1.) / 2.;
            const auto pole = 1. / (std::complex<double>(0, wn) - omega[j]);
            const std::complex<double> expected = bose ? (i == 0 ? -1. : omega[j]*pole)
                : std::string(type) == "anomalous" ? -omega[j]*pole : pole;
            EXPECT_NEAR(k.K(2*i,j), expected.real(), 1.e-14);
            EXPECT_NEAR(k.K(2*i+1,j), expected.imag(), 1.e-14);
        }
    }
}
TEST(MaxEntKernel, BosonicTimeIsEvenAndHasThermalZeroLimit) {
    auto p = parameters("bosonic", "time");
    kernel_probe k(p, {-3., 0., 3., 1.e-15});
    for (int i = 0; i < 4; ++i) {
        const double tau = 2.*i/3;
        const double expected = 1.5 * std::cosh(3*(1-tau))/std::sinh(3.);
        EXPECT_NEAR(k.K(i,0), expected, 1.e-14);
        EXPECT_NEAR(k.K(i,2), expected, 1.e-14);
        EXPECT_DOUBLE_EQ(k.K(i,1), .5); EXPECT_NEAR(k.K(i,3), .5, 1.e-14);
    }
    kernel_probe large(p, {-1000., 1000.});
    EXPECT_TRUE(std::isfinite(large.K(0,0))); EXPECT_EQ(large.K(0,0), large.K(0,1));
}
TEST(MaxEntKernel, BosonicSymmetricFrequencyHasContinuousZeroLimit) {
    auto p = parameters("bosonic", "frequency"); p["PARTICLE_HOLE_SYMMETRY"] = true;
    kernel_probe k(p, {0., 1.e-12});
    EXPECT_DOUBLE_EQ(k.K(0,0), -1.); EXPECT_DOUBLE_EQ(k.K(0,1), -1.);
    for (int i = 1; i < 4; ++i) EXPECT_EQ(k.K(i,0), 0.);
}
TEST(MaxEntKernel, OddComplexDataFailBeforeWritingPastTheMatrix) {
    const auto p = parameters("fermionic", "frequency", 5);
    EXPECT_THROW((kernel_probe(p, {1., 2.})), std::invalid_argument);
}
