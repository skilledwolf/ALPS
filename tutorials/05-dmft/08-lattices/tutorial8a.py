# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2012-2013 by Jakub Imriska  <jimriska@phys.ethz.ch>
# 
# ALPS Project: https://alps.comp-phys.org/
# SPDX-License-Identifier: MIT
# 
# ****************************************************************************

import pyalps
import matplotlib.pyplot as plt
import pyalps.plot
from pyalps.run_io import execute, write_run_file


#prepare the input parameters and run the simulation
for u in [3.]:
  for b in [6.]:
    name = 'hybrid_DOS_beta_'+str(b)+'_U_'+str(u)
    run = write_run_file(name+'.toml', overwrite=True,
        parameters={
            'BETA' : b,          # inverse temperature
            'MU' : 0.0,          # chemical potential corresponding to half-filling
            'U' : u,             # Hubbard interaction
            'FLAVORS' : 2,       # corresponds to spin up/down
            'SITES' : 1,         # number of sites in the impurity
            'H' : 0.0,           # there is no magnetic field
            'H_INIT' : 0.05,     # we set initial field to split spin up/down in order to trigger AF phase
            'ANTIFERROMAGNET' : True,   # allow Neel order
            'SYMMETRIZATION' : False,   # do not enforce paramagnetic solution
            'NMATSUBARA' : 500,         # number of Matsubara frequencies
            'N' : 500,                  # bins in imaginary time
            'CONVERGED' : 0.005,        # criterion for convergency
            'SWEEPS' : 10000,           # max. number of sweeps in a single iteration
            'THERMALIZATION' : 500,     # number of thermalization sweeps
            'N_MEAS' : 5000,            # number of Monte Carlo steps between measurements
            'N_HISTOGRAM_ORDERS' : 50,  # histogram size
            'EPS_0' : 0.,               # potential shift for the flavor 0
            'EPS_1' : 0.,               # potential shift for the flavor 1
            'EPSSQ_0' : 4.,             # the second moment of the bandstructure for the flavor 0
            'EPSSQ_1' : 4.,             # the second moment of the bandstructure for the flavor 1
        },
        input={'dos': 'DOS/DOS_Square_GRID4000'},  # file with the density of states
        output={'results': name+'.h5'},
        # The self-consistency runs in Matsubara frequencies and starts from the local
        # non-interacting Green's function. The hybridization solver receives Delta.
        execution={'solver': 'hybridization',
                   'max_iterations': 10,  # max. number of self-consistency iterations
                   'time_limit': 60,      # max. time spent in solver in a single iteration in seconds
                   'seed': 0})
    execute('dmft', run)

listobs=['0']  # we look only at flavor=0
    
data = pyalps.loadMeasurements(pyalps.getResultFiles(pattern='hybrid_DOS*h5'), respath='/simulation/results/G_tau', what=listobs, verbose=True)
for d in pyalps.flatten(data):
    d.x = d.x*d.props["BETA"]/float(d.props["N"])
    d.props['label'] = r'$\beta=$'+str(d.props['BETA'])
plt.figure()

plt.xlabel(r'$\tau$')
plt.ylabel(r'$G_{flavor=0}(\tau)$')
plt.title('DMFT-08, DOS-based approach: Hubbard model on the square lattice')
pyalps.plot.plot(data)
plt.legend()
plt.show()
