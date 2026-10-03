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
from pyalps.run_io import execute, write_run_file

#prepare the input parameters
# For more precise calculations we propose to you to:
#   enhance the time_limit and max_iterations, and lower CONVERGED
run = write_run_file('parm_hyb.toml', overwrite=True,
    parameters={
        'ANTIFERROMAGNET'     : False,
        'CONVERGED'           : 0.0025,
        'FLAVORS'             : 2,
        'H'                   : 0.,
        'H_INIT'              : 0.0,
        'MU'                  : 0.,
        'N'                   : 1000,
        'NMATSUBARA'          : 1000,
        'N_MEAS'              : 10000,
        'N_HISTOGRAM_ORDERS'  : 50,
        'SITES'               : 1,
        'SYMMETRIZATION'      : True,
        'U'                   : 3.,
        't'                   : 0.707106781186547,
        'SWEEPS'              : 2500,
        'THERMALIZATION'      : 500,
        'BETA'                : 32.
    },
    output={'results': 'parm_hyb.h5'},
    execution={'solver': 'hybridization', 'max_iterations': 12, 'time_limit': 600, 'seed': 0})

#run the simulation
execute('dmft', run)
