// SPDX-License-Identifier: MIT
#include "ising.hpp"
#include <gtest/gtest.h>
#include <cmath>

TEST(IsingTutorial, SamplesTheMeasuredFerromagneticHamiltonian)
{
    alps::params p;
    p["L"] = 8;
    p["T"] = 2.;
    p["SEED"] = 42;
    p["THERMALIZATION"] = 1000;
    p["SWEEPS"] = 20000;
    ising_sim simulation(p);
    for (int i = 0; i < 21000; ++i) {
        simulation.update();
        simulation.measure();
    }
    const auto results = simulation.collect_results();
    // Differentiate Z = (2 cosh(beta))^L + (2 sinh(beta))^L.
    const double t = std::tanh(0.5);
    const double expected = -(t + std::pow(t, 7)) / (1. + std::pow(t, 8));
    EXPECT_EQ(results["Energy"].count(), 20000u);
    // Fixed seed, broad tolerance: this detects the wrong ensemble, not
    // convergence of the error estimator or a specific random trajectory.
    EXPECT_NEAR(results["Energy"].mean<double>(), expected, 0.06);
}
