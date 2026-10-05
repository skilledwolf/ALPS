// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
#ifdef ALPS_NDIM_HEISENBERG
#include "o_n_model/ndim_spin.hpp"
using simulation_type = ndim_spin_sim<4>;
#else
#include "heisenberg.hpp"
using simulation_type = heisenberg_sim;
#endif
#include <alps/hdf5/archive.hpp>
#include <filesystem>
#include <stdexcept>

int main(int argc, char* argv[]) {
    alps::params parameters;
    parameters["L"] = 3;
    parameters["T"] = 2.;
    parameters["THERMALIZATION"] = 11;
    parameters["SWEEPS"] = 701;
    parameters["SEED"] = 42;
#if defined(ALPS_NDIM_HEISENBERG) || defined(ALPS_LATTICE_HEISENBERG)
    parameters["LATTICE"] = "simple cubic lattice";
    if (argc == 2) parameters["LATTICE_LIBRARY"] = std::string(argv[1]);
#endif
    auto different = parameters;
    different["T"] = 7.;
    different["THERMALIZATION"] = 3;
    different["SWEEPS"] = 99;
    simulation_type uninterrupted(parameters), stopped(parameters), resumed(different);
    for (int i = 0; i < 712; ++i) {
        uninterrupted.update(); uninterrupted.measure();
        if (i < 148) { stopped.update(); stopped.measure(); }
    }
    auto filename = std::string(parameters.exists("LATTICE") ? "lattice-heisenberg-" : "chain-heisenberg-")
                  + std::to_string(stopped.get_measurements().at("Magnetization")->size()) + ".h5";
    stopped.save(filename);
    resumed.load(filename);
    for (int i = 148; i < 712; ++i) { resumed.update(); resumed.measure(); }
    auto expected = uninterrupted.collect_results();
    auto check_results = [&] {
        auto actual = resumed.collect_results();
        for (auto const& entry : expected) {
            auto const& result = actual.at(entry.first);
            if (entry.second != result || entry.second.store().count() != result.store().count()
                    || entry.second.stderror() != result.stderror())
                throw std::runtime_error("resuming changed " + entry.first);
        }
    };
    check_results();
    if (uninterrupted.fraction_completed() != resumed.fraction_completed()
            || uninterrupted.get_random()() != resumed.get_random()())
        throw std::runtime_error("resuming changed Heisenberg progress or random stream");
    // Reject bad application state before changing any restored base state.
    {
        alps::hdf5::archive archive(filename, "a");
        archive["/parameters"] << different;
        archive["/simulation/realizations/0/clones/0/checkpoint/sweeps"] << -1;
    }
    bool rejected = false;
    try { resumed.load(filename); } catch (std::invalid_argument const&) { rejected = true; }
    std::filesystem::remove(filename);
    if (!rejected || resumed.get_parameters()["T"].as<double>() != 2.
            || resumed.fraction_completed() != 1.)
        throw std::runtime_error("invalid Heisenberg checkpoint changed application state");
    check_results();
}
