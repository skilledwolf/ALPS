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
# For more precise calculations we propose to enhance the SWEEPS
for u,j in [[1.8,0.45],[2.2,0.55],[2.8,0.7]]:
    name = 'parm_u_'+str(u)+'_j_'+str(j)
    run = write_run_file(name+'.toml', overwrite=True,
        parameters={
            'CONVERGED'           : 0.001,
            'FLAVORS'             : 4,
            'H'                   : 0.,
            'H_INIT'              : 0.,
            'MU'                  : 0.,
            'N'                   : 500,
            'NMATSUBARA'          : 500,
            'N_MEAS'              : 2000,
            'N_HISTOGRAM_ORDERS'  : 50,
            'SYMMETRIZATION'      : True,
            'SWEEPS'              : 10000,
            'BETA'                : 30.,
            'THERMALIZATION'      : 500,
            'U'                   : u,
            'J'                   : j,
            't0'                  : 0.5,
            't1'                  : 1.
        },
        output={'results': name+'.h5'},
        execution={'solver': 'hybridization', 'max_iterations': 15, 'time_limit': 600, 'seed': 0})
    execute('dmft', run)

listobs = ['0', '2']   # flavor 0 is SYMMETRIZED with 1, flavor 2 is SYMMETRIZED with 3
    
data = pyalps.loadMeasurements(pyalps.getResultFiles(pattern='parm_u_*h5'), respath='/simulation/results/G_tau', what=listobs, verbose=True)
for d in pyalps.flatten(data):
    d.x = d.x*d.props["BETA"]/float(d.props["N"])
    d.y = -d.y
    d.props['label'] = r'$U=$'+str(d.props['U'])+'; flavor='+str(d.props['observable'][len(d.props['observable'])-1])
plt.figure()
plt.yscale('log')
plt.xlabel(r'$\tau$')
plt.ylabel(r'$G_{flavor}(\tau)$')
plt.title('DMFT-05: Orbitally Selective Mott Transition on the Bethe lattice')
pyalps.plot.plot(data)
plt.legend()
plt.show()
