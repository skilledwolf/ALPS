# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2009-2010 by Matthias Troyer <troyer@phys.ethz.ch> 
#                            Jan Gukelberger
#                            Brigitte Surer
# 
# ALPS Project: https://alps.comp-phys.org/
# SPDX-License-Identifier: MIT
# 
# ****************************************************************************

from solution.ising import Simulation as NativeSimulation, main


class Simulation(NativeSimulation):
    # Implement only the physics; parameters, RNG, statistics and output are
    # shared with the solution so the exercise uses the same native interfaces.
    def step(self):
        # Attempt L*L random-site Metropolis proposals with periodic neighbors.
        # Propose keeping or flipping the spin with equal probability.
        raise NotImplementedError('Implement the Metropolis sweep')

    def observables(self):
        # Return the per-site [E, m, abs(m), m**2, m**4] sample.
        raise NotImplementedError('Implement the observables')


if __name__ == '__main__':
    main(simulation=Simulation)
