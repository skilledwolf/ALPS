# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2010 by Jan Gukelberger <gukelberger@phys.ethz.ch> 
# 
# ALPS Project: https://alps.comp-phys.org/
# SPDX-License-Identifier: MIT
# 
# ****************************************************************************

import pyalps
import numpy as np
import matplotlib.pyplot as plt
import pyalps.plot

#prepare the input parameters
parms= []
for m in [20,40,60]:
    parms.append({ 
        'LATTICE'                   : 'open chain lattice with special edges 32',
        'MODEL'                     : 'spin',
        'local_S0'                  : '0.5',
        'local_S1'                  : '1',
        'CONSERVED_QUANTUMNUMBERS'  : 'N,Sz',
        'Sz_total'                  : 0,
        'J'                         : 1,
        'SWEEPS'                    : 4,
        'NUMBER_EIGENVALUES'        : 1,
        'MAXSTATES'                 : m
       })

#write the input file and run the simulation
from pyalps.run_io import execute, write_run_files
runs = [dict(parameters=p, input=dict(lattice_library='my_lattices.xml'),
             output=dict(results=f'parm_spin_one_multiple.task{i+1}.out.h5'))
        for i, p in enumerate(parms)]
files = execute('dmrg', write_run_files('parm_spin_one_multiple', runs, overwrite=True))

#load all measurements for all states
data = pyalps.loadEigenstateMeasurements(files)

# print properties of the eigenvector for each run:
for run in data:
    for s in run:
        print(s.props['observable'], ' : ', s.y[0])
