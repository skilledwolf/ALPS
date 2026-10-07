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
parms = [ { 
        'LATTICE'                   : 'open chain lattice with special edges',
        'MODEL'                     : 'spin',
        'local_S0'                  : '0.5',
        'local_S1'                  : '1',
        'CONSERVED_QUANTUMNUMBERS'  : 'N,Sz',
        'Sz_total'                  : 0,
        'J'                         : 1,
        'SWEEPS'                    : 4,
        'NUMBER_EIGENVALUES'        : 1,
        'MAXSTATES'                 : 100
       } ]

#write the input file and run the simulation
from pyalps.run_io import execute, write_run_files
runs = [dict(parameters=p, input=dict(lattice_library='my_lattice.xml'),
             output=dict(results=f'parm_spin_one.task{i+1}.out.h5'))
        for i, p in enumerate(parms)]
files = execute('dmrg', write_run_files('parm_spin_one', runs, overwrite=True))

#load all measurements for all states
data = pyalps.loadEigenstateMeasurements(files)

# print properties of the eigenvector:
for s in data[0]:
    print(s.props['observable'], ' : ', s.y[0])

# load and plot iteration history
iter = pyalps.loadMeasurements(files,
                               what=['Iteration Energy','Iteration Truncation Error'])

plt.figure()
pyalps.plot.plot(iter[0][0])
plt.title('Iteration history of ground state energy (S=1)')
plt.ylabel('$E_0$')
plt.xlabel('iteration')

plt.figure()
pyalps.plot.plot(iter[0][1])
plt.title('Iteration history of truncation error (S=1)')
plt.yscale('log')
plt.ylabel('error')
plt.xlabel('iteration')

plt.show()
