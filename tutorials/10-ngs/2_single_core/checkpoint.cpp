// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
#include "ising.hpp"
#include <cmath>
#include <filesystem>
#include <stdexcept>

int main() {
    alps::params parameters;
    parameters["L"] = 7;
    parameters["T"] = 2.;
    parameters["THERMALIZATION"] = 11;
    parameters["SWEEPS"] = 701;
    parameters["SEED"] = 42;
    ising_sim uninterrupted(parameters), stopped(parameters), resumed(parameters);
    for (int i = 0; i < 712; ++i) {
        uninterrupted.update(); uninterrupted.measure();
        if (i < 148) { stopped.update(); stopped.measure(); }
    }
    auto filename = "mcbase-ising-continuation.h5";
    stopped.save(filename); // Restore an unfinished native ALEA batch.
    resumed.load(filename);
    for (int i = 148; i < 712; ++i) { resumed.update(); resumed.measure(); }
    auto expected = uninterrupted.collect_results();
    double u = std::tanh(1. / double(parameters["T"]));
    double energy = -(u + std::pow(u, int(parameters["L"]) - 1))
                  / (1. + std::pow(u, int(parameters["L"])));
    auto const& measured = expected.at("Energy");
    if (!(std::abs(measured.mean()(0) - energy) <= 6. * measured.stderror()(0)))
        throw std::runtime_error("Ising energy disagrees with the Boltzmann distribution");
    auto check_results = [&] {
        auto actual = resumed.collect_results();
        for (auto const& entry : expected) {
            auto const& restored = actual.at(entry.first);
            if (entry.second != restored || entry.second.store().count() != restored.store().count()
                    || entry.second.stderror() != restored.stderror())
                throw std::runtime_error("resuming changed " + entry.first);
        }
    };
    check_results();
    if (uninterrupted.fraction_completed() != resumed.fraction_completed())
        throw std::runtime_error("resuming changed simulation progress");
    if (uninterrupted.get_random()() != resumed.get_random()())
        throw std::runtime_error("resuming changed the random stream");
    {
        alps::hdf5::archive archive(filename, "a");
        auto different = parameters;
        different["T"] = 7.;
        archive["/parameters"] << different;
        archive["/simulation/realizations/0/clones/0/checkpoint/sweeps"] << -1;
    }
    bool rejected = false;
    try { resumed.load(filename); } catch (std::invalid_argument const&) { rejected = true; }
    std::filesystem::remove(filename);
    if (!rejected || resumed.get_parameters()["T"].as<double>() != 2.
            || resumed.fraction_completed() != 1.)
        throw std::runtime_error("invalid Ising checkpoint changed application state");
    check_results();

}
