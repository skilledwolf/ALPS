# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2010 by Brigitte Surer <surerb@phys.ethz.ch> 
#               2012 by Jakub Imriska  <jimriska@phys.ethz.ch>
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
run = write_run_file('parm_int.toml', overwrite=True,
    parameters={
        'ANTIFERROMAGNET'         : False,
        'CONVERGED'               : 0.0025,
        'FLAVORS'                 : 2,
        'H'                       : 0.,
        'H_INIT'                  : 0.,
        'MU'                      : 0.,
        'N'                       : 500,
        'NMATSUBARA'              : 500,
        'NMATSUBARA_MEASUREMENTS' : 18,
        'NSELF'                   : 5000,
        'MEASUREMENT_PERIOD'      : 10,
        'SITES'                   : 1,
        'SYMMETRIZATION'          : True,
        'U'                       : 3.,
        't'                       : 0.707106781186547,
        'RECALC_PERIOD'           : 3000,
        'SWEEPS'                  : 100000000,
        'THERMALIZATION'          : 1000,
        'ALPHA'                   : -0.01,
        'HISTOGRAM_MEASUREMENT'   : True,
        'BETA'                    : 32.
    },
    output={'results': 'parm_int.h5'},
    execution={'solver': 'interaction', 'max_iterations': 12, 'time_limit': 120, 'seed': 0})

#run the simulation
execute('dmft', run)
