// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
#include "heisenberg.hpp"
#include <cmath>
#include <stdexcept>

int main() {
    alps::params parameters;
    parameters["L"] = 2;
    parameters["T"] = 2.;
    parameters["THERMALIZATION"] = 1024;
    parameters["SWEEPS"] = 32768;
    parameters["SEED"] = 42;
    heisenberg_sim simulation(parameters);
    simulation.run([] { return false; });
    auto energy = simulation.collect_results().at("Energy");
    // Two periodic bonds give Z proportional to sinh(2 beta)/(2 beta).
    double coupling = 2. / double(parameters["T"]);
    double exact = 1. / coupling - 1. / std::tanh(coupling);
    if (!(std::abs(energy.mean()(0) - exact) <= 6. * energy.stderror()(0)))
        throw std::runtime_error("Heisenberg Metropolis samples the wrong coupling");
}
