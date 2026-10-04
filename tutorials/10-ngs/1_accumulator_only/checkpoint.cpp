// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
#include "ising.hpp"
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
    for (auto const& entry : expected) {
        auto const& restored = actual.at(entry.first);
        if (entry.second != restored || entry.second.store().count() != restored.store().count()
                || entry.second.stderror() != restored.stderror())
            throw std::runtime_error("resuming changed " + entry.first);
    }
    if (uninterrupted.fraction_completed() != resumed.fraction_completed())
        throw std::runtime_error("resuming changed simulation progress");
}
