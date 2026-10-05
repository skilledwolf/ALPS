/*****************************************************************************
 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations
 *
 * ALPS Libraries
 *
 * Copyright (C) 2003 by Brigitte Surer
 *                       and Jan Gukelberger
 *
 * ALPS Project: https://alps.comp-phys.org/
 * SPDX-License-Identifier: MIT
 *
 *****************************************************************************/

#include "simulation.hpp"

// Implement the physics here. The solution uses the same parameters, RNG,
// native ALEA statistics and HDF5 output defined in simulation.hpp.
void Simulation::step() {
    // Attempt L*L random-site Metropolis flips with periodic neighbors.
    ...
}

std::vector<double> Simulation::observables() const {
    // Return the per-site {E, m, abs(m), m*m, m*m*m*m} sample.
    ...
}

int main(int argc, char** argv) { return run_tutorial(argc, argv); }
