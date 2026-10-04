// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#include "hirschfyesim.h"
#include "run_config.h"
#include <alps/alea/hdf5.hpp>
#include <alps/random/buffered_rng.h>
#include <cmath>
#include <filesystem>
#include <iostream>

namespace {
void require(bool value, char const* message) {
  if (!value) throw std::runtime_error(message);
}
void close(double actual, double expected, char const* message, double tolerance = 1.e-12) {
  require(std::abs(actual - expected) < tolerance, message);
}
template<class F> void rejects(F operation, char const* message) {
  try { operation(); } catch (std::exception const&) { return; }
  throw std::runtime_error(message);
}
alps::run_configuration configuration(double u = 0., int thermalization = 2, int samples = 3) {
  alps::run_configuration run;
  run.parameters["BETA"] = 2.; run.parameters["U"] = u;
  run.parameters["N"] = 4; run.parameters["NMATSUBARA"] = 4;
  run.parameters["EPSSQ_0"] = 0.; run.parameters["EPSSQ_1"] = 0.;
  run.parameters["THERMALIZATION"] = thermalization; run.parameters["SWEEPS"] = samples;
  run.parameters = alps::dmft::prepare_hirschfye_parameters(run.parameters);
  run.execution["seed"] = 19; run.execution["bins"] = 16;
  return run;
}
matsubara_green_function_t free_green() {
  matsubara_green_function_t green(4, 1, 2);
  for (int flavor = 0; flavor < 2; ++flavor)
    for (int i = 0; i < 4; ++i) green(i, flavor) = {0., -2. / ((2 * i + 1) * std::acos(-1.))};
  return green;
}
struct fixture : HirschFyeRun {
  using HirschFyeRun::HirschFyeRun;
  using HirschFyeRun::spins;
  int configuration_sign() const { return sign; }
  void check_wide_progress() {
    thermalization_sweeps = 0; total_sweeps = std::uint64_t{1} << 60;
    sweeps = total_sweeps - 1;
    require(work_done() < 1., "floating-point progress ended a wide run before its final sweep");
    ++sweeps;
    require(work_done() == 1., "wide run did not finish at its requested sweep count");
  }
  double random() { return random_01(); }
  auto raw() const { return measurements[1].result(); }
  void sample(double sample_sign, std::array<double, 4> const& values) {
    measurements[0] << sample_sign;
    alps::alea::column<double> green(6);
    for (int i = 0; i < 4; ++i) green(i) = sample_sign * values[i];
    green(4) = sample_sign * (-1. - values[0]); green(5) = sample_sign;
    measurements[1] << green; measurements[2] << green;
  }
  itime_green_function_t exact_configuration(std::vector<int> const& configuration) const {
    dense_matrix up, down;
    const auto weight_sign = update_from_zero(up, Green0_up, configuration, lambda)
                           * update_from_zero(down, Green0_down, configuration, -lambda);
    require(weight_sign == 1, "physical atomic configuration has negative weight");
    itime_green_function_t tau(N + 1, 1, 2);
    green_vector_from_matrix(tau, up, down);
    return tau;
  }
  void check_rebuild_bounds() {
    // U=0 updates leave the injected off-diagonal error untouched.
    N_check = 1; check_counter = 0; Green_up(0, N - 1) += 1.;
    dostep();
    require(N_check == 1, "precision deterioration reduced rebuild interval to zero");
    close(norm_max(Green_up), .5, "rebuild did not restore the free matrix");
    N_check = std::numeric_limits<std::uint64_t>::max() / 2 + 1;
    check_counter = N_check - 1;
    dostep();
    require(N_check == std::numeric_limits<std::uint64_t>::max(), "rebuild interval overflowed");
  }
};
struct draws {
  double site, acceptance;
  bool first = true;
  double operator()() { if (first) { first = false; return site; } return acceptance; }
};

void numerical_contracts() {
  close(hirschfye_lambda(2., 1., 4), std::acosh(std::exp(.25)), "HS coupling differs from exact formula");
  close(hirschfye_lambda(2., 1.e-20, 4), std::sqrt(5.e-21), "small HS coupling lost significance", 1.e-24);
  rejects([] { hirschfye_lambda(2., 1.e300, 4); }, "overflowing HS exponent was accepted");
  dense_matrix matrix(2, 2); matrix.clear(); matrix(0, 1) = -7.;
  close(norm_max(matrix), 7., "matrix norm omitted a non-first-column element");
  matrix(1, 1) = std::numeric_limits<double>::quiet_NaN();
  require(std::isinf(norm_max(matrix)), "matrix norm hid a nonfinite element");

  dense_matrix bare(1, 1), restored(1, 1); restored(0, 0) = 42.;
  bare(0, 0) = 2.;
  rejects([&] { update_from_zero(restored, bare, {1}, std::log(2.)); }, "singular LAPACK solve was accepted");
  close(restored(0, 0), 42., "failed solve modified its destination");
  bare(0, 0) = 3.;
  require(update_from_zero(restored, bare, {1}, std::log(2.)) == -1, "negative determinant sign was lost");
  close(restored(0, 0), -3., "negative determinant solve is incorrect");
  // A=[[0,1],[1,0]] has negative determinant solely from the pivot swap.
  bare.resize(2, 2); bare(0, 0) = 2.; bare(1, 1) = 2.; bare(0, 1) = -1.; bare(1, 0) = -1.;
  require(update_from_zero(restored, bare, {1, 1}, std::log(2.)) == -1, "LU permutation sign was lost");

  dense_matrix up(1, 1), down(1, 1); up(0, 0) = 1.; down(0, 0) = 2.;
  std::vector<int> spin{1}; int sign = 1;
  draws reject{0., .8};
  update_single_spin(reject, up, down, spin, std::log(2.), sign);
  require(spin[0] == 1 && sign == 1, "negative ratio -2 used a heat-bath probability above 2/3");
  draws accept{0., .5};
  update_single_spin(accept, up, down, spin, std::log(2.), sign);
  require(spin[0] == -1 && sign == -1, "accepted negative ratio did not flip configuration sign");
  close(down(0, 0), -1., "negative-ratio rank-one update is incorrect");
  draws reverse_reject{0., .5};
  update_single_spin(reverse_reject, up, down, spin, std::log(2.), sign);
  require(spin[0] == -1 && sign == -1, "reverse negative ratio -1/2 used probability above 1/3");
  draws reverse_accept{0., .25};
  update_single_spin(reverse_accept, up, down, spin, std::log(2.), sign);
  require(spin[0] == 1 && sign == 1, "reverse accepted update lost sign");
  close(down(0, 0), 2., "reverse rank-one update did not recover the original matrix");
}

void exact_atomic_science() {
  auto run = configuration(1.);
  fixture simulation(run, free_green());
  const double coupling = hirschfye_lambda(2., 1., 4);
  std::array<long double, 4> mean{};
  long double partition = 0.;
  for (int bits = 0; bits < 16; ++bits) {
    std::vector<int> spins(4);
    int sum = 0;
    for (int i = 0; i < 4; ++i) { spins[i] = (bits & (1 << i)) ? 1 : -1; sum += spins[i]; }
    // Independent closed-form determinant product for the atomic bath.
    const long double weight = (1. + std::cosh(coupling * sum)) / 2.;
    auto tau = simulation.exact_configuration(spins);
    for (int i = 0; i < 4; ++i) mean[i] -= weight * tau(i, 0);
    partition += weight;
  }
  for (int time = 0; time < 4; ++time) {
    // Four-state trace of U(n_up-1/2)(n_down-1/2), independent of HF matrices.
    const double exact = -std::cosh((.5 * time - 1.) / 2.) / (2. * std::cosh(.5));
    close(mean[time] / partition, exact, "finite-U spin enumeration differs from exact atomic trace");
  }
  // A rank-one accepted update must equal the independently rebuilt matrix.
  dense_matrix bare(4, 4), up, down;
  for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) bare(i, j) = i >= j ? .5 : -.5;
  std::vector<int> spins{1, -1, 1, -1};
  update_from_zero(up, bare, spins, coupling); update_from_zero(down, bare, spins, -coupling);
  int sign = 1;
  const auto old_up = up, old_down = down;
  const double ratio = (1. + std::cosh(2. * coupling)) / 2.;
  draws boundary{.3, ratio / (1. + ratio) + 1.e-8};
  update_single_spin(boundary, up, down, spins, coupling, sign);
  require(spins[1] == -1, "atomic heat-bath rejection differs from determinant-weight oracle");
  require(norm_max(dense_matrix(up - old_up)) == 0. && norm_max(dense_matrix(down - old_down)) == 0., "rejected update changed matrices");
  draws accepted{.3, 0.};
  update_single_spin(accepted, up, down, spins, coupling, sign);
  dense_matrix exact_up, exact_down;
  update_from_zero(exact_up, bare, spins, coupling); update_from_zero(exact_down, bare, spins, -coupling);
  close(norm_max(dense_matrix(up - exact_up)), 0., "rank-one up update differs from exact rebuild");
  close(norm_max(dense_matrix(down - exact_down)), 0., "rank-one down update differs from exact rebuild");
  require(sign == 1, "positive atomic weight changed sign");
}

void streams_and_counts() {
  auto run = configuration();
  for (int rank : {0, 3}) {
    fixture simulation(run, free_green(), rank);
    alps::buffered_rng<boost::mt19937> engine;
    engine.seed(19 + rank);
    boost::variate_generator<alps::buffered_rng<boost::mt19937>&, boost::uniform_real<>> old(engine, boost::uniform_real<>{});
    for (auto spin : simulation.spins)
      require(spin == (old() < .5 ? 1 : -1), "initial auxiliary spins changed the legacy seeded stream");
    for (int i = 0; i < 11000; ++i) close(simulation.random(), old(), "uniform stream differs across legacy buffer boundaries", 1.e-16);
  }
  for (int warmup : {0, 2}) {
    run = configuration(0., warmup, 3);
    fixture simulation(run, free_green());
    for (int step = 0; step < warmup; ++step) simulation.dostep();
    require(simulation.collect_results().at("Sign").count() == 0, "warmup contributed measurements");
    for (int step = 0; step < 3; ++step) simulation.dostep();
    close(simulation.work_done(), 1., "work progress includes warmup or wall time");
    auto results = simulation.collect_results();
    require(results.at("Sign").count() == 3 && results.at("G_meas_up").count() == 3, "measurement count differs from requested completed sweeps");
    const auto green = HirschFyeRun::get_result(results, run.parameters);
    for (int spin = 0; spin < 2; ++spin) for (int time = 0; time <= 4; ++time) {
      close(green.second(time, spin), -.5, "free physical time Green function is incorrect");
      close(green.second.error(time, spin), 0., "free Green function has nonzero statistical error");
    }
    for (int i = 0; i < 4; ++i) require(std::abs(green.first(i, 0) - free_green()(i, 0)) < 1.e-12, "free frequency Green function is incorrect");
  }
  fixture bounds(configuration(), free_green()); bounds.check_rebuild_bounds();
  bounds.check_wide_progress();
  // Finite input validation does not imply causal/positive-weight G0. With
  // internal G0(0)=3 in both channels exactly one determinant is negative,
  // independently of the randomly selected initial auxiliary spin.
  run = configuration(std::log(1.25)); run.parameters["N"] = 1;
  auto unphysical = free_green();
  for (int spin = 0; spin < 2; ++spin) unphysical(0, spin).real(-2.5);
  fixture negative(run, unphysical);
  require(negative.configuration_sign() == -1, "initialization assumed positive weight for every accepted finite bath");
}

void signed_statistics() {
  auto run = configuration();
  fixture simulation(run, free_green());
  auto empty = simulation.collect_results().at("G_meas_up");
  require(empty.size() == 5 && empty.count() == 0, "empty signed Green result lost its shape");
  fixture single(run, free_green()); single.sample(-1., {-.3, -.4, -.4, -.4});
  auto once = single.collect_results().at("G_meas_up");
  close(once.mean()(0), -.3, "single-bin signed Green mean was discarded");
  require(once.stderror()(0) == std::numeric_limits<double>::infinity(), "single-bin signed error implies independent samples");
  std::array<long double, 6> sum{};
  constexpr int samples = 173;
  for (int i = 0; i < samples; ++i) {
    const double sign = i % 7 ? 1. : -1.;
    const std::array<double, 4> value = {-.4 - .08 * std::sin(i), -.3 - .01 * std::cos(i), -.2 - .03 * std::sin(2 * i), -.3 - .02 * std::cos(2 * i)};
    simulation.sample(sign, value);
    for (int j = 0; j < 4; ++j) sum[j] += sign * value[j];
    sum[4] += sign * (-1. - value[0]); sum[5] += sign;
  }
  const auto raw = simulation.raw();
  const auto results = simulation.collect_results();
  auto const& result = results.at("G_meas_up");
  require(raw == simulation.raw() && result.count() == samples && result.store().count() == raw.store().count(), "signed projection consumed state or lost partial-bin weights");
  long double count2 = 0.;
  for (auto weight : raw.store().count()) count2 += static_cast<long double>(weight) * weight;
  std::array<long double, 5> mean{};
  std::vector<std::array<long double, 5>> pseudo(raw.num_batches());
  for (std::size_t bin = 0; bin < raw.num_batches(); ++bin) {
    const auto weight = raw.store().count()(bin);
    if (!weight) continue;
    for (int j = 0; j < 5; ++j) {
      const auto leave = (sum[j] - raw.store().batch()(j, bin)) / (sum[5] - raw.store().batch()(5, bin));
      pseudo[bin][j] = (samples * sum[j] / sum[5] - (samples - weight) * leave) / weight;
      mean[j] += weight * pseudo[bin][j] / samples;
    }
  }
  const auto covariance = (result.cov() * result.count2() / (double(samples) * samples)).eval();
  for (int i = 0; i < 5; ++i) {
    close(result.mean()(i), mean[i], "signed Green mean differs from independent pseudovalues");
    for (int j = 0; j < 5; ++j) {
      long double squared = 0.;
      for (std::size_t bin = 0; bin < pseudo.size(); ++bin)
        squared += raw.store().count()(bin) * (pseudo[bin][i] - mean[i]) * (pseudo[bin][j] - mean[j]);
      const auto expected = squared / (samples - count2 / samples) * count2 / (samples * samples);
      close(covariance(i, j), expected, "signed endpoint covariance differs from independent bin oracle");
    }
    close(result.stderror()(i) * result.stderror()(i), covariance(i, i), "signed error disagrees with covariance diagonal");
    close(covariance(i, 4), -covariance(i, 0), "endpoint closure lost anticorrelated errors");
  }
  close(result.mean()(0) + result.mean()(4), -1., "signed physical endpoint closure failed");
  const auto green = HirschFyeRun::get_result(results, run.parameters);
  for (int time = 0; time <= 4; ++time) close(green.second.error(time, 0), result.stderror()(time), "derived Green function lost native errors");
  fixture zero(run, free_green()); zero.sample(1., {-.3, -.4, -.4, -.4}); zero.sample(-1., {-.3, -.4, -.4, -.4});
  const auto before = zero.raw();
  rejects([&] { zero.collect_results(); }, "zero pooled sign was accepted");
  require(before == zero.raw(), "failed ratio analysis changed live measurements");
  const auto file = std::filesystem::path("hirschfye-statistics.h5");
  require(!std::filesystem::exists(file), "test output already exists");
  {
    alps::hdf5::archive archive(file.string(), "w");
    alps::alea::hdf5_serializer serializer(archive, "/simulation/results");
    serialize(serializer, "G_meas_up", result);
  }
  {
    alps::hdf5::archive archive(file.string(), "r");
    alps::alea::hdf5_serializer serializer(archive, "/simulation/results");
    alps::alea::batch_result<double> restored;
    deserialize(serializer, "G_meas_up", restored);
    require(restored == result, "native result codec changed physical endpoint bins");
  }
  std::filesystem::remove(file);
}
}
int main() {
  try {
    numerical_contracts(); exact_atomic_science(); streams_and_counts(); signed_statistics();
    return 0;
  } catch (std::exception const& error) { std::cerr << error.what() << '\n'; return 1; }
}
