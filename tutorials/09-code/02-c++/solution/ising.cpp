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

#include "../simulation.hpp"

void Simulation::step() {
    for (size_t s=0; s<size_t(L)*L; ++s) {
        const int i = randint(L), j = randint(L);
        if (random()<.5) continue; // symmetric no-change proposal prevents parity trapping
        const int neighbors = spins[(i+L-1)%L][j] + spins[(i+1)%L][j]
                            + spins[i][(j+L-1)%L] + spins[i][(j+1)%L];
        const int e = -spins[i][j]*neighbors;
        // Only uphill moves need an exponential, avoiding overflow at large beta.
        if (e > 0 || random() < std::exp(beta*(2*e))) spins[i][j] = -spins[i][j];
    }
}

std::vector<double> Simulation::observables() const {
    double energy = 0, magnetization = 0;
    for (int i=0; i<L; ++i) for (int j=0; j<L; ++j) {
        energy -= spins[i][j]*(spins[(i+1)%L][j] + spins[i][(j+1)%L]);
        magnetization += spins[i][j];
    }
    const double m = magnetization/(double(L)*L);
    return {energy/(double(L)*L), m, std::abs(m), m*m, m*m*m*m};
}

int main(int argc, char** argv) { return run_tutorial(argc, argv); }
