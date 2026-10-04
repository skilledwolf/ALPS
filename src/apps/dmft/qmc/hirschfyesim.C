/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2009 by Emanuel Gull <gull@phys.columbia.edu>
 *                              Philipp Werner <werner@itp.phys.ethz.ch>,
 *                              Sebastian Fuchs <fuchs@theorie.physik.uni-goettingen.de>
 *                              Matthias Troyer <troyer@comp-phys.org>
 *
 *
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
 *
 *****************************************************************************/

#include "hirschfyesim.h"
#include "fouriertransform.h"
#include <alps/alea/transform.hpp>
#include <alps/random/seed.h>
#include <boost/numeric/bindings/lapack/driver/gesv.hpp>
#include <exception>

// Include every matrix element: column-major begin1() alone visits only one
// column and can miss a loss of precision elsewhere.
double norm_max(dense_matrix const& matrix) {
  double maximum = 0.;
  for (std::size_t column = 0; column < matrix.size2(); ++column)
    for (std::size_t row = 0; row < matrix.size1(); ++row) {
      const double value = std::abs(matrix(row, column));
      if (!std::isfinite(value)) return std::numeric_limits<double>::infinity();
      maximum = std::max(maximum, value);
    }
  return maximum;
}

int update_from_zero(dense_matrix& green, dense_matrix const& green0,
                     std::vector<int> const& spins, double lambda) {
  const int n = spins.size();
  dense_matrix matrix(n, n), restored(green0);
  for (int column = 0; column < n; ++column) {
    const double exponential = std::exp(lambda * spins[column]);
    for (int row = 0; row < n; ++row)
      matrix(row, column) = -green0(row, column) * (exponential - 1.);
    matrix(column, column) += exponential;
  }
  boost::numeric::ublas::vector<fortran_int_t> pivots(n);
  const auto status = boost::numeric::bindings::lapack::gesv(matrix, pivots, restored);
  if (status != 0 || !std::isfinite(norm_max(restored)))
    throw std::runtime_error("Hirsch-Fye Green-function matrix solve failed");
  int sign = 1;
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(matrix(i, i)))
      throw std::runtime_error("Hirsch-Fye Green-function matrix solve failed");
    if (matrix(i, i) < 0.) sign *= -1;
    if (pivots(i) != i + 1) sign *= -1;
  }
  green.swap(restored);
  return sign;
}

HirschFyeRun::HirschFyeRun(alps::run_configuration const& run,
                         matsubara_green_function_t const& g0, int rank)
  : parameters(run.parameters),
    thermalization_sweeps(parameters["THERMALIZATION"].as<std::uint64_t>()),
    total_sweeps(parameters["SWEEPS"].as<std::uint64_t>()),
    beta(parameters["BETA"].as<double>()), N(parameters["N"].as<int>()),
    lambda(hirschfye_lambda(beta, parameters["U"].as<double>(), N)),
    tolerance(parameters["TOLERANCE"].as<double>()),
    Green0_up(N, N), Green_up(N, N), Green0_down(N, N), Green_down(N, N),
    spins(N), bare_green_tau(N + 1, 1, 2), green_tau(N + 1, 1, 2) {
  if (parameters["SITES"].as<int>() != 1 || g0.nsite() != 1 || g0.nflavor() != 2 ||
      g0.nfreq() != parameters["NMATSUBARA"].as<unsigned int>())
    throw std::invalid_argument("Hirsch-Fye requires a one-site, two-flavor Green function of NMATSUBARA values");
  alps::seed_with_sequence(random_01.engine(),
    run.execution["seed"].as<std::uint32_t>() + static_cast<std::uint32_t>(rank));
  boost::shared_ptr<FourierTransformer> fourier;
  FourierTransformer::generate_transformer(parameters, fourier);
  fourier->backward_ft(bare_green_tau, g0);
  for (int spin = 0; spin < 2; ++spin)
    for (int time = 0; time <= N; ++time)
      bare_green_tau(time, spin) *= -1.; // Internal matrix convention is positive.
  green_matrix_from_vector(bare_green_tau, Green0_up, Green0_down);
  for (auto& spin : spins) spin = random_01() < .5 ? 1 : -1;
  sign = update_from_zero(Green_up, Green0_up, spins, lambda)
       * update_from_zero(Green_down, Green0_down, spins, -lambda);
  const auto bins = run.execution["bins"].as<std::size_t>();
  measurements = {alps::alea::batch_acc<double>(1, bins, 1),
                  alps::alea::batch_acc<double>(N + 2, bins, 1),
                  alps::alea::batch_acc<double>(N + 2, bins, 1)};
}

double HirschFyeRun::work_done() const {
  const auto completed = sweeps > thermalization_sweeps ? sweeps - thermalization_sweeps : 0;
  return completed >= total_sweeps ? 1.
    : std::min(completed / double(total_sweeps), std::nextafter(1., 0.));
}

void HirschFyeRun::dostep() {
  ++sweeps;
  ++check_counter;
  for (int i = 0; i < N; ++i)
    update_single_spin(random_01, Green_up, Green_down, spins, lambda, sign);
  if (sweeps > thermalization_sweeps) {
    // The scheduler discarded its inclusive warm-up-transition measurement.
    // Recording strictly afterwards preserves exactly SWEEPS samples.
    green_vector_from_matrix(green_tau, Green_up, Green_down);
    measurements[0] << double(sign);
    alps::alea::column<double> sample(N + 2);
    for (int spin = 0; spin < 2; ++spin) {
      for (int time = 0; time < N; ++time)
        sample(time) = -sign * green_tau(time, spin);
      // Affine endpoint closure belongs in each sample, before statistics.
      sample(N) = sign * (green_tau(0, spin) - 1.);
      sample(N + 1) = sign;
      measurements[spin + 1] << sample;
    }
  }
  if (check_counter >= N_check) {
    dense_matrix restored_up, restored_down;
    const int restored_sign = update_from_zero(restored_up, Green0_up, spins, lambda)
                            * update_from_zero(restored_down, Green0_down, spins, -lambda);
    dense_matrix difference_up(Green_up), difference_down(Green_down);
    difference_up -= restored_up;
    difference_down -= restored_down;
    Green_up.swap(restored_up);
    Green_down.swap(restored_down);
    sign = restored_sign;
    if (norm_max(difference_up) > tolerance || norm_max(difference_down) > tolerance)
      N_check = std::max<std::uint64_t>(1, N_check / 2);
    else
      N_check = N_check > std::numeric_limits<std::uint64_t>::max() / 2
        ? std::numeric_limits<std::uint64_t>::max() : N_check * 2;
    check_counter = 0;
  }
}

namespace {
struct green_ratio final : alps::alea::transformer<double> {
  explicit green_ratio(std::size_t size) : size(size) {}
  std::size_t in_size() const override { return size + 1; }
  std::size_t out_size() const override { return size; }
  alps::alea::column<double> operator()(alps::alea::column<double> const& value) const override {
    if (value(size) == 0. || !std::isfinite(value(size)))
      throw std::domain_error("Hirsch-Fye signed estimate has zero or nonfinite average sign");
    return value.head(size) / value(size);
  }
  std::size_t size;
};
constexpr std::array<char const*, 3> names = {"Sign", "G_meas_up", "G_meas_down"};
}

HirschFyeRun::results_type HirschFyeRun::collect_results(alps::alea::reducer const* reduction) const {
  results_type results;
  // The fixed three-entry registry cannot differ in names or cadence. Reduce
  // every raw block before any recipient-only ratio can throw.
  for (std::size_t i = 0; i < measurements.size(); ++i) {
    auto result = measurements[i].result();
    if (reduction) result.reduce(*reduction);
    if (result.valid()) results.emplace(names[i], std::move(result));
  }
  std::exception_ptr failure;
  try {
    for (std::size_t i = 1; i < names.size() && !results.empty(); ++i) {
      auto& joint = results.at(names[i]);
      green_ratio ratio(joint.size() - 1);
      const auto occupied = (joint.store().count().array() != 0).count();
      if (occupied > 1) {
        joint = alps::alea::transform(alps::alea::jackknife_prop{}, ratio, joint);
      } else {
        alps::alea::batch_data<double> values(ratio.out_size(), joint.num_batches());
        values.count() = joint.store().count();
        if (joint.count()) {
          const auto mean = ratio(joint.mean());
          for (std::size_t bin = 0; bin < joint.num_batches(); ++bin)
            if (values.count()(bin)) values.batch().col(bin) = values.count()(bin) * mean;
        }
        joint = alps::alea::batch_result<double>(values);
      }
    }
  } catch (...) { failure = std::current_exception(); }
  if (reduction && reduction->get_max(bool(failure))) {
    if (failure) std::rethrow_exception(failure);
    throw std::runtime_error("Hirsch-Fye result analysis failed on the recipient");
  }
  if (failure) std::rethrow_exception(failure);
  return results;
}

std::pair<matsubara_green_function_t, itime_green_function_t>
HirschFyeRun::get_result(results_type const& results, alps::params const& parameters) {
  const auto slices = parameters["N"].as<unsigned int>();
  itime_green_function_t tau(slices + 1, 1, 2);
  matsubara_green_function_t omega(parameters["NMATSUBARA"].as<unsigned int>(), 1, 2);
  std::vector<double> densities(2);
  for (int spin = 0; spin < 2; ++spin) {
    auto const& result = results.at(names[spin + 1]);
    if (!result.count() || result.size() != slices + 1)
      throw std::invalid_argument("Hirsch-Fye Green-function result has no samples or incorrect shape");
    const auto mean = result.mean();
    const auto error = result.stderror();
    for (unsigned int time = 0; time <= slices; ++time) {
      tau(time, spin) = mean(time);
      tau.error(time, spin) = error(time);
    }
    densities[spin] = -tau(slices, spin);
  }
  require_finite(tau, "Hirsch-Fye G(tau)");
  boost::shared_ptr<FourierTransformer> fourier;
  FourierTransformer::generate_transformer_U(parameters, fourier, densities);
  fourier->forward_ft(tau, omega);
  require_finite(omega, "Hirsch-Fye G(iw)");
  return {omega, tau};
}

void HirschFyeRun::green_matrix_from_vector(itime_green_function_t const& bare,
                                            dense_matrix& up, dense_matrix& down) const {
  for (int column = 0; column < N; ++column) {
    for (int row = column; row < N; ++row) {
      up(row, column) = bare(row - column, 0);
      down(row, column) = bare(row - column, 1);
    }
    for (int row = column + 1; row < N; ++row) {
      up(column, row) = -bare(N - (row - column), 0);
      down(column, row) = -bare(N - (row - column), 1);
    }
  }
}

void HirschFyeRun::green_vector_from_matrix(itime_green_function_t& tau,
                                            dense_matrix const& up, dense_matrix const& down) const {
  tau.clear();
  for (int column = 0; column < N; ++column) {
    for (int row = column; row < N; ++row) {
      tau(row - column, 0) += up(row, column);
      tau(row - column, 1) += down(row, column);
      if (row > column) {
        tau(N - (row - column), 0) -= up(column, row);
        tau(N - (row - column), 1) -= down(column, row);
      }
    }
  }
  for (int time = 0; time < N; ++time) {
    tau(time, 0) /= N;
    tau(time, 1) /= N;
  }
}
