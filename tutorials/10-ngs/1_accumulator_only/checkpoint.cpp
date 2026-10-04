// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
#include "ising.hpp"
#include <cmath>
#include <filesystem>
#include <stdexcept>

int main() {
    auto filename = "ising-continuation.h5";
    alps::params parameters;
    parameters["L"] = 7;
    parameters["T"] = 2.;
    parameters["THERMALIZATION"] = 11;
    parameters["SWEEPS"] = 701;
    parameters["SEED"] = 42;
    ising_sim uninterrupted(parameters), stopped(parameters), resumed(parameters);
    for (int i = 0; i < 712; ++i) {
        uninterrupted.update();
        uninterrupted.measure();
        if (i < 148) { stopped.update(); stopped.measure(); }
    }
    stopped.save(filename); // Deliberately inside an unfinished ALEA batch.
    resumed.load(filename);
    std::filesystem::remove(filename);
    for (int i = 148; i < 712; ++i) { resumed.update(); resumed.measure(); }
    auto expected = uninterrupted.collect_results(), actual = resumed.collect_results();
    // Exact finite periodic-chain energy from the Ising transfer matrix.
    double u = std::tanh(1. / double(parameters["T"]));
    double energy = -(u + std::pow(u, int(parameters["L"]) - 1))
                  / (1. + std::pow(u, int(parameters["L"])));
    auto const& measured = expected.at("Energy");
    if (!(std::abs(measured.mean()(0) - energy) <= 6. * measured.stderror()(0)))
        throw std::runtime_error("Ising energy disagrees with the Boltzmann distribution");
    for (auto const& entry : expected) {
        auto const& restored = actual.at(entry.first);
        if (entry.second != restored || entry.second.store().count() != restored.store().count()
                || entry.second.stderror() != restored.stderror())
            throw std::runtime_error("resuming changed " + entry.first);
    }
    if (uninterrupted.fraction_completed() != resumed.fraction_completed())
        throw std::runtime_error("resuming changed simulation progress");
}
