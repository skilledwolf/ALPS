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
# For more precise simulation we propose to you to:
#   lower CONVERGED (to 0.0003) and TOLERANCE (to 0.0001) and increase the time_limit
for b in [6., 12.]:
    run = write_run_file('parm_beta_'+str(b)+'.toml', overwrite=True,
        parameters={
            'ANTIFERROMAGNET'     : True,
            'CONVERGED'           : 0.005,
            'FLAVORS'             : 2,
            'H'                   : 0.,
            'H_INIT'              : 0.05,
            'MU'                  : 0.,
            'N'                   : 16,
            'NMATSUBARA'          : 500,
            'SITES'               : 1,
            'SYMMETRIZATION'      : False,
            'TOLERANCE'           : 0.001,
            'U'                   : 3.,
            't'                   : 0.707106781186547,
            'SWEEPS'              : 1000000,
            'THERMALIZATION'      : 10000,
            'BETA'                : b
        },
        output={'results': 'parm_beta_'+str(b)+'.h5'},
        execution={'solver': 'hirschfye', 'max_iterations': 10,
                   'time_limit': 20, 'seed': 0})
    execute('dmft', run)

listobs=['0', '1']
    
data = pyalps.loadMeasurements(pyalps.getResultFiles(pattern='parm_beta_*h5'), respath='/simulation/results/G_tau', what=listobs)
for d in pyalps.flatten(data):
    d.x = d.x*d.props["BETA"]/float(d.props["N"])
    d.props['label'] = r'$\beta=$'+str(d.props['BETA'])+'; flavor='+str(d.props['observable'][len(d.props['observable'])-1])
plt.figure()
plt.xlabel(r'$\tau$')
plt.ylabel(r'$G_{flavor}(\tau)$')
plt.title('DMFT-07: Neel transition for the Hubbard model on the Bethe lattice\n(using the Hirsch-Fye impurity solver)')
pyalps.plot.plot(data)
plt.legend()
plt.show()
