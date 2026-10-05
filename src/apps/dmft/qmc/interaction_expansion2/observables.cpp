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

#include "interaction_expansion.hpp"
#include <alps/alea/transform.hpp>
#include <exception>

namespace {
struct signed_ratio final : alps::alea::transformer<double> {
  explicit signed_ratio(std::size_t components) : components(components) {}
  std::size_t in_size() const override { return components + 1; }
  std::size_t out_size() const override { return components; }
  alps::alea::column<double> operator()(alps::alea::column<double> const& value) const override {
    if (value(components) == 0.)
      throw std::domain_error("CT-INT signed estimate has zero average sign");
    return value.head(components) / value(components);
  }
  std::size_t components;
};
}

void InteractionExpansion::initialize_observables()
{
  measurements.clear();
  auto add = [&](std::string const& name, std::size_t components, bool signed_value = false) {
    measurements.emplace(name, measurement{signed_value,
      alps::alea::batch_acc<double>(components + signed_value, num_bins, 1)});
  };
  add("Sign", 1);
  add("PertOrder", n_flavors);
  for (unsigned int flavor = 0; flavor < n_flavors; ++flavor) {
    if (measurement_method == selfenergy_measurement_itime_rs) {
      for (unsigned int i = 0; i < n_site; ++i) {
        for (unsigned int j = 0; j < n_site; ++j)
          add("W_" + std::to_string(flavor) + "_" + std::to_string(i) + "_" + std::to_string(j), n_self + 1, true);
        auto name = "density_" + std::to_string(flavor);
        if (n_site > 1) name += "_" + std::to_string(i);
        add(name, 1, true);
      }
    } else {
      for (unsigned int k = 0; k < n_site; ++k) {
        auto suffix = std::to_string(flavor) + "_" + std::to_string(k) + "_" + std::to_string(k);
        add("Wk_real_" + suffix, n_matsubara_measurements, true);
        add("Wk_imag_" + suffix, n_matsubara_measurements, true);
      }
      add("densities_" + std::to_string(flavor), n_site, true);
    }
  }
  if (measurement_method == selfenergy_measurement_itime_rs && n_flavors == 2) {
    for (unsigned int i = 0; i < n_site; ++i) {
      auto suffix = std::to_string(i);
      add("Sz_" + suffix, 1, true);
      add("Sz2_" + suffix, 1, true);
      add("Sz0_Sz" + suffix, 1, true);
    }
  }
  add("densities", n_flavors, true);
  add("density_correlation", 1, true);
  add("n_i n_j", std::size_t(n_flavors) * n_flavors, true);
  add("VertexInsertion", 1);
  add("VertexRemoval", 1);
}

void InteractionExpansion::record_measurement(std::string const& name, std::valarray<double> const& value)
{
  auto& entry = measurements.at(name);
  alps::alea::column<double> sample(value.size() + entry.signed_value);
  if (sample.size() != entry.accumulator.size()) throw alps::alea::size_mismatch();
  for (std::size_t i = 0; i < value.size(); ++i) sample(i) = value[i];
  // The numerator is already sign-weighted by the physical estimator. Its
  // denominator is recorded at this same cadence, never from another stream.
  if (entry.signed_value) sample(value.size()) = sign;
  entry.accumulator << sample;
}

void InteractionExpansion::record_measurement(std::string const& name, double value)
{
  record_measurement(name, std::valarray<double>(value, 1));
}

InteractionExpansion::results_type InteractionExpansion::collect_results(alps::alea::reducer const* reduction) const
{
  if (reduction) {
    auto agree = [&](std::int64_t value) {
      auto maximum = reduction->get_max(value);
      auto minimum = -reduction->get_max(-value);
      if (maximum != minimum)
        throw std::runtime_error("CT-INT measurement registries differ between replicas");
    };
    agree(measurements.size());
    for (auto const& entry : measurements) {
      agree(entry.first.size());
      for (unsigned char byte : entry.first) agree(byte);
      agree(entry.second.signed_value);
    }
  }
  results_type results;
  // Every raw block is reduced before any recipient-only ratio analysis.
  // Independent replicas must contribute their raw numerators and signs.
  for (auto const& entry : measurements) {
    auto result = entry.second.accumulator.result();
    if (reduction) result.reduce(*reduction);
    if (result.valid()) results.emplace(entry.first, std::move(result));
  }
  std::exception_ptr failure;
  try {
    for (auto& entry : results) {
      if (!measurements.at(entry.first).signed_value) continue;
      auto const& joint = entry.second;
      signed_ratio ratio(joint.size() - 1);
      auto occupied = (joint.store().count().array() != 0).count();
      if (occupied > 1) {
        entry.second = alps::alea::transform(alps::alea::jackknife_prop{}, ratio, joint);
      } else {
        // Empty results preserve shape/count without evaluating a ratio.
        // One occupied bin has a mean, but no independent error estimate.
        alps::alea::batch_data<double> values(ratio.out_size(), joint.num_batches());
        values.count() = joint.store().count();
        if (joint.count()) {
          auto mean = ratio(joint.mean());
          for (std::size_t i = 0; i < joint.num_batches(); ++i)
            if (values.count()(i)) values.batch().col(i) = values.count()(i) * mean;
        }
        entry.second = alps::alea::batch_result<double>(values);
      }
    }
  } catch (...) { failure = std::current_exception(); }
  if (reduction && reduction->get_max(bool(failure))) {
    if (failure) std::rethrow_exception(failure);
    throw std::runtime_error("CT-INT result analysis failed on the recipient");
  }
  if (failure) std::rethrow_exception(failure);
  return results;
}

void InteractionExpansion::measure_observables()
{
  record_measurement("Sign", sign);
  if (measurement_method == selfenergy_measurement_matsubara)
    compute_W_matsubara();
  else if (measurement_method == selfenergy_measurement_itime_rs)
    compute_W_itime();
  measure_densities();
  std::valarray<double> pert_order(n_flavors);
  for (unsigned int i = 0; i < n_flavors; ++i) {
    assert(num_rows(M[i].matrix()) == num_cols(M[i].matrix()));
    pert_order[i] = num_rows(M[i].matrix());
  }
  record_measurement("PertOrder", pert_order);
}
