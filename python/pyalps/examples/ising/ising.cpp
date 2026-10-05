// Copyright (C) 2010-2012 by Lukas Gamper
//               2026      by the ALPS collaboration
// SPDX-License-Identifier: MIT

#include "ising.hpp"

#include <alps/hdf5/archive.hpp>
#include <alps/ngs.hpp>

#include <algorithm>
#include <stdexcept>

ising_sim::ising_sim(parameters_type const & parameters,
                     std::size_t seed_offset)
    : alps::mcbase(parameters, seed_offset),
      total_sweeps_(parameters.exists("SWEEPS") ? parameters["SWEEPS"].as<std::size_t>() : 10) {
    if (!total_sweeps_) throw std::invalid_argument("SWEEPS must be positive");
    measurements.emplace("Magnetization", std::make_shared<alps::alea::batch_acc<double>>(1, 64));
}

void ising_sim::update() {
    state_ = random() < 0.5 ? -1.0 : 1.0;
    ++sweeps_;
}

void ising_sim::measure() {
    *measurement("Magnetization") << alps::alea::make_adapter(state_);
}

double ising_sim::fraction_completed() const {
    return std::min(1.0, static_cast<double>(sweeps_) / total_sweeps_);
}

void ising_sim::save(alps::hdf5::archive & archive) const {
    alps::mcbase::save(archive);
    archive["checkpoint/sweeps"] << sweeps_;
    archive["checkpoint/state"] << state_;
}

void ising_sim::load(alps::hdf5::archive & archive) {
    parameters_type restored_parameters;
    archive["/parameters"] >> restored_parameters;
    const auto total = restored_parameters.value_or<std::size_t>("SWEEPS", 10);
    std::size_t sweeps;
    double state;
    archive["checkpoint/sweeps"] >> sweeps;
    archive["checkpoint/state"] >> state;
    std::vector<std::uint64_t> counts;
    archive["measurements/Magnetization/batch/count"] >> counts;
    std::uint64_t count = 0;
    for (auto value : counts) {
        if (value > sweeps - count) throw std::invalid_argument("invalid measurement count");
        count += value;
    }
    if (!total || sweeps > total || (state != -1 && state != 1) || count != sweeps)
        throw std::invalid_argument("invalid simulation checkpoint");
    alps::mcbase::load(archive);
    total_sweeps_ = total;
    sweeps_ = sweeps;
    state_ = state;
}
