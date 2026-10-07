# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2009-2010 by Matthias Troyer <troyer@phys.ethz.ch> 
# 
# ALPS Project: https://alps.comp-phys.org/
# SPDX-License-Identifier: MIT
# 
# ****************************************************************************

import pyalps
import matplotlib.pyplot as plt
import pyalps.plot
import numpy as np

#prepare the input parameters
parms = [{ 
          'LATTICE'                   : "chain lattice", 
          'MODEL'                     : "spin",
          'CONSERVED_QUANTUMNUMBERS'  : 'Sz',
          'local_S'                   : 1,
          'J'                         : 1,
          'L'                         : 8
        }]

#write the input file and run the simulation
from pyalps.run_io import execute, write_run_files
runs = [dict(parameters=p, output=dict(results=f'parm6a.task{i+1}.out.h5'))
        for i, p in enumerate(parms)]
files = execute('fulldiag', write_run_files('parm6a', runs, overwrite=True))

#run the evaluation and load all the plots
data = pyalps.evaluateFulldiagVersusT(files,DELTA_T=0.1, T_MIN=0.1, T_MAX=10.0)

#make plot
for s in pyalps.flatten(data):
  plt.figure()
  plt.title("Antiferromagnetic Heisenberg chain")
  pyalps.plot.plot(s)
plt.show()

