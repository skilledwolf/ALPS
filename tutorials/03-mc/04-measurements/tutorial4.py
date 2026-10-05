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

#prepare the input parameters
parms = [{ 
          'LATTICE'                   : "square lattice", 
          'MODEL'                     : "spin",
          'MEASURE[Correlations]'     : True,
          'MEASURE[Structure Factor]' : True,
          'MEASURE[Green Function]'   : True,
          'local_S'                   : 0.5,
          'T'                         : 0.3,
          'J'                         : 1 ,
          'THERMALIZATION'            : 10000,
          'SWEEPS'                    : 500000,
          'L'                         : 4,
          'h'                         : 0.1
        }]

#write the input file and run the simulation
from pyalps.run_io import execute, write_run_files
runs = [dict(parameters=p, output=dict(results=f'parm4.task{i+1}.out.h5'))
        for i, p in enumerate(parms)]
files = execute('dirloop_sse', write_run_files('parm4', runs, baseseed=42, overwrite=True))

#load the magnetization and collect it as function of field h
data = pyalps.loadMeasurements(files)

# print all measurements
for s in pyalps.flatten(data):
  if len(s.x)==1:
    print(s.props['observable'], ' : ', s.y[0])
  else:
    for (x,y) in zip(s.x,s.y):
      print(s.props['observable'], x, ' : ', y)
