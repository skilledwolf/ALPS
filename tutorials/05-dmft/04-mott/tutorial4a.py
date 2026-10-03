# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2010      by Brigitte Surer <surerb@phys.ethz.ch> 
#               2012-2013 by Jakub Imriska  <jimriska@phys.ethz.ch>
# 
# ALPS Project: https://alps.comp-phys.org/
# SPDX-License-Identifier: MIT
# 
# ****************************************************************************

import pyalps
import numpy as np
import matplotlib.pyplot as plt
import pyalps.plot
from pyalps.run_io import execute, write_run_file


#prepare the input parameters and run the simulation
# For more precise calculations change the parameters
#   enlarge SWEEPS,
#   CONVERGED (to 0.0002-0.001),
#   raise N and NMATSUBARA (to 1000)
for u in [4.,5.,6.,8.]:
    run = write_run_file('parm_u_'+str(u)+'.toml', overwrite=True,
        parameters={
            'ANTIFERROMAGNET'         : False,
            'CONVERGED'               : 0.001,
            'FLAVORS'                 : 2,
            'H'                       : 0.,
            'H_INIT'                  : 0.,
            'MU'                      : 0.,
            'N'                       : 500,
            'NMATSUBARA'              : 500,
            'N_MEAS'                  : 1000,
            'N_HISTOGRAM_ORDERS'      : 50,
            'SITES'                   : 1,
            'SYMMETRIZATION'          : True,
            't'                       : 1.,
            'SWEEPS'                  : int(1500*u),
            'BETA'                    : 20.0,
            'THERMALIZATION'          : 500,
            'U'                       : u
        },
        output={'results': 'parm_u_'+str(u)+'.h5'},
        execution={'solver': 'hybridization', 'max_iterations': 20, 'time_limit': 600, 'seed': 0})
    execute('dmft', run)

listobs=['0']   # we look at only one flavor, as they are SYMMETRIZED
    
data = pyalps.loadMeasurements(pyalps.getResultFiles(pattern='parm_u_*h5'), respath='/simulation/results/G_tau', what=listobs, verbose=True)

for d in pyalps.flatten(data):
    d.x = d.x*d.props["BETA"]/float(d.props["N"])
    d.y = -d.y
    d.props['label'] = r'$U=$'+str(d.props['U'])
plt.figure()
plt.yscale('log')
plt.xlabel(r'$\tau$')
plt.ylabel(r'$G_{flavor=0}(\tau)$')
plt.title('DMFT-04: Mott-insulator transition for the Hubbard model on the Bethe lattice')
pyalps.plot.plot(data)
plt.legend()
plt.show()
