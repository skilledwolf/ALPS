# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2025 by ALPS Collaboration
#
# ALPS Project: https://alps.comp-phys.org/
# SPDX-License-Identifier: MIT
# 
# ****************************************************************************

import pyalps
import numpy as np

#prepare the input parameters
parms = []
sz_tot = [1,2]
for sz in sz_tot:
    parms.append( {
        'LATTICE'                   : "open chain lattice", 
        'MODEL'                     : "spin",
        'local_S'                   : '1',
        'CONSERVED_QUANTUMNUMBERS'  : 'N,Sz',
        'Sz_total'                  : sz,
        'J'                         : 1,
        'SWEEPS'                    : 5,
        'L'                         : 64,
        'MAXSTATES'                 : 300,
        'NUMBER_EIGENVALUES'        : 1
       } )

#write the input file and run the simulation
from pyalps.run_io import execute, write_run_files
runs = [dict(parameters=p, output=dict(results=f'parm_spin_one_triplet.task{i+1}.out.h5'))
        for i, p in enumerate(parms)]
files = execute('dmrg', write_run_files('parm_spin_one_triplet', runs, overwrite=True))

#load all measurements for all states
data = pyalps.loadEigenstateMeasurements(files)

# print results:
energies = {}
for run in data:
    print('S_z =', run[0].props['Sz_total'])
    for s in run:
        print('\t', s.props['observable'], ':', s.y[0])
        if s.props['observable'] == 'Energy':
            sz = s.props['Sz_total']
            energies[sz] = s.y[0]

print('Gap:', energies[sz_tot[1]]-energies[sz_tot[0]])
