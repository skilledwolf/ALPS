// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
// Exact configuration oracles for the retained CT-HYB estimators. The signed
// zero-frequency correction being tested adapts ALPS reconciliation commit
// 90a36c0d4c64ee7308dcb0651264197ffdf36449, hybridization/hyblocal.cpp.
// Expected values come from integrating piecewise constant occupations and
// tracing H = U*n_up*n_down - MU*(n_up+n_down) over its four occupation states.
// No Markov-chain convergence or postprocessed ratio uncertainty is claimed.

#include <gtest/gtest.h>
#include "hybridization/hyblocal.hpp"
#include <cmath>
#include <utility>
#include <vector>

namespace {
alps::params parameters(double beta, double interaction = 0., double mu = 0.) {
  alps::params p;
  p["BETA"] = beta;
  p["N_ORBITALS"] = 2;
  p["N_TAU"] = 16;
  p["U"] = interaction;
  p["MU"] = mu;
  return p;
}

// Real part of integral_start^end exp(i*omega*tau) d tau, with the zero
// mode treated as the interval length. Unlike the production implementation,
// this oracle uses the scalar sine antiderivative and explicitly splits wraps.
double interval_integral(double start, double end, double omega) {
  return omega == 0. ? end - start
                    : (std::sin(omega * end) - std::sin(omega * start)) / omega;
}

double occupation_integral(const std::vector<std::pair<double, double>>& intervals,
                           double beta, std::size_t mode) {
  const double omega = 2. * std::acos(-1.) * mode / beta;
  double value = 0.;
  for (const auto& interval : intervals) {
    if (interval.second >= interval.first)
      value += interval_integral(interval.first, interval.second, omega);
    else {
      value += interval_integral(interval.first, beta, omega);
      value += interval_integral(0., interval.second, omega);
    }
  }
  return value;
}
}

TEST(CthybDensityOracle, EmptyAndFullLinesHaveSignedZeroModeOnly) {
  for (double beta : {2., 5.}) {
    for (double sign : {-1., 1.}) {
      for (bool occupied : {false, true}) {
        for (std::size_t modes : {std::size_t(1), std::size_t(6)}) {
          SCOPED_TRACE(beta);
          SCOPED_TRACE(sign);
          SCOPED_TRACE(occupied);
          SCOPED_TRACE(modes);
          local_configuration configuration(parameters(beta), 0);
          configuration.set_zero_order_orbital_occupied(0, occupied);
          std::vector<double> density(2, 0.);
          configuration.measure_density(density, sign);
          EXPECT_DOUBLE_EQ(density[0], occupied ? sign : 0.);
          EXPECT_DOUBLE_EQ(configuration.measure_nn(0, 0), occupied ? 1. : 0.);
          std::vector<double> frequency(modes, 0.);
          configuration.measure_nnw(0, frequency, sign);
          EXPECT_DOUBLE_EQ(frequency[0], occupied ? sign * beta : 0.);
          for (std::size_t mode = 1; mode < modes; ++mode)
            EXPECT_DOUBLE_EQ(frequency[mode], 0.);
        }
      }
    }
  }
}

TEST(CthybDensityOracle, SegmentModesMatchExactOccupationIntegrals) {
  for (double beta : {2., 5.}) {
    for (bool winding : {false, true}) {
      SCOPED_TRACE(beta);
      SCOPED_TRACE(winding);
      const std::vector<std::pair<double, double>> intervals = winding
        ? std::vector<std::pair<double, double>>{{.8 * beta, .2 * beta},
                                                 {.3 * beta, .5 * beta}}
        : std::vector<std::pair<double, double>>{{.1 * beta, .4 * beta},
                                                 {.6 * beta, .85 * beta}};
      local_configuration configuration(parameters(beta), 0);
      for (const auto& interval : intervals)
        configuration.insert_segment(segment(interval.first, interval.second), 0);
      EXPECT_NEAR(configuration.segment_density(0), winding ? .6 : .55, 1.e-12);
      for (double sign : {-1., 1.}) {
        SCOPED_TRACE(sign);
        // A nonzero initial value also checks the estimator's additive contract.
        std::vector<double> frequency(6, .125);
        configuration.measure_nnw(0, frequency, sign);
        for (std::size_t mode = 0; mode < frequency.size(); ++mode)
          EXPECT_NEAR(frequency[mode],
                      .125 + sign * occupation_integral(intervals, beta, mode),
                      1.e-12 * beta);
      }
    }
  }
}

TEST(CthybDensityOracle, AtomicDensityCorrelationsMatchTheFourStateTrace) {
  for (double beta : {2., 5.}) {
    for (const auto& model : {std::pair<double, double>{0., 0.},
                             std::pair<double, double>{1.7, .4}}) {
      const double interaction = model.first;
      const double mu = model.second;
      const double single_weight = std::exp(beta * mu);
      const double double_weight = std::exp(beta * (2. * mu - interaction));
      const double partition = 1. + 2. * single_weight + double_weight;
      const double occupation = (single_weight + double_weight) / partition;
      const double double_occupation = double_weight / partition;
      SCOPED_TRACE(beta);
      SCOPED_TRACE(interaction);
      SCOPED_TRACE(mu);
      for (double sign : {-1., 1.}) {
        SCOPED_TRACE(sign);
        std::vector<double> density(2, 0.);
        std::vector<double> diagonal(6, 0.);
        std::vector<double> off_diagonal(6, 0.);
        double equal_time_pair = 0.;
        for (int state = 0; state < 4; ++state) {
          const bool up = state & 1;
          const bool down = state & 2;
          const double energy = interaction * up * down - mu * (up + down);
          const double probability = std::exp(-beta * energy) / partition;
          local_configuration configuration(parameters(beta, interaction, mu), 0);
          configuration.set_zero_order_orbital_occupied(0, up);
          configuration.set_zero_order_orbital_occupied(1, down);
          std::vector<double> state_density(2, 0.);
          configuration.measure_density(state_density, sign);
          for (std::size_t orbital = 0; orbital < density.size(); ++orbital)
            density[orbital] += probability * state_density[orbital];
          equal_time_pair += probability * configuration.measure_nn(0, 1);
          // Match production measure_nnw's n_j(0) occupation gate. For an
          // atomic state n_i(tau) is constant and all nonzero modes vanish.
          for (int second_orbital = 0; second_orbital < 2; ++second_orbital) {
            if (configuration.density(second_orbital, 0.) == 0.) continue;
            std::vector<double> state_frequency(6, 0.);
            configuration.measure_nnw(0, state_frequency, sign);
            auto& destination = second_orbital == 0 ? diagonal : off_diagonal;
            for (std::size_t mode = 0; mode < destination.size(); ++mode)
              destination[mode] += probability * state_frequency[mode];
          }
        }
        for (double value : density)
          EXPECT_NEAR(value / sign, occupation, 1.e-12);
        EXPECT_NEAR(equal_time_pair, double_occupation, 1.e-12);
        EXPECT_NEAR(diagonal[0] / sign, beta * occupation, 1.e-12);
        EXPECT_NEAR(off_diagonal[0] / sign, beta * double_occupation, 1.e-12);
        for (std::size_t mode = 1; mode < diagonal.size(); ++mode) {
          EXPECT_DOUBLE_EQ(diagonal[mode], 0.);
          EXPECT_DOUBLE_EQ(off_diagonal[mode], 0.);
        }
      }
    }
  }
}
