 #############################################################################/
 #
 # ALPS Project: Algorithms and Libraries for Physics Simulations
 #
 # ALPS Libraries
 #
 # Copyright (C) 2012 by Hartmut Hafermann <hafermann@cpht.polytechnique.fr>
 #
 #
 # ALPS Project: https://alps.comp-phys.org/
 # SPDX-License-Identifier: MIT
 #
 #############################################################################/

 # This tutorial is a minimal example illustrating the use of the python interface
 # to hybridization expansion solver.
 #
 # The hybridization function is chosen such that the impurity model is equivalent
 # to a correlated site couple to a single bath site with coupling V=1 and energy
 # epsilon=0 (Delta(tau)=-V**2/2=const.) . It can therefore be compared to the
 # exact result obtained by exact diagonalization (see subdirectory ED).
 #
 # Run this script as:
 # alpspython tutorial1.py
 #

import pyalps.cthyb as cthyb # the solver module

# specify solver parameters
parms={
'SWEEPS'              : 100000000,
'THERMALIZATION'      : 1000,
'N_MEAS'              : 50,
'N_HISTOGRAM_ORDERS'  : 50,
'N_ORBITALS'          : 2,
'U'                   : 4.0,
'MU'                  : 2.0,
'N_TAU'               : 1000,
'BETA'                : 45.,
}

# Write a simple (constant) hybridization function to file.
# This corresponds to the impurity site with on-site interaction U
# coupled to a single noninteracting bath site with hybridization
# V=1 at energy epsilon=0: Delta(tau)=-V^2/2=const.

f=open("delta.dat","w")
for i in range(parms["N_TAU"]+1):
    f.write("%i %f %f\n"%(i,-0.5,-0.5))
f.close()

# solve the impurity model for 60 seconds; results are written to results.out.h5
# and, in human-readable form, to text files in the current directory
run=cthyb.prepare(parms, input={'delta': "delta.dat"},
                  output={'results': "results.out.h5", 'text': True},
                  execution={'time_limit': 60, 'seed': 0})
cthyb.solve(run)




