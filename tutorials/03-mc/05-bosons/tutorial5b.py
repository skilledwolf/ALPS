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

#prepare the input parameters
parms = []
for t in [0.045, 0.05, 0.0525, 0.055, 0.0575, 0.06, 0.065]:
  for l in [4, 6, 8]:
    parms.append(
        { 
          'LATTICE'        : "square lattice", 
          'MODEL'          : "boson Hubbard",
          'T'              : 0.05,
          'L'              : l ,
          't'              : t ,
          'U'              : 1.0 ,
          'mu'             : 0.5, 
          'NONLOCAL'       : False ,
          'Nmax'           : 2 ,
          'THERMALIZATION' : 15000,
          'SWEEPS'         : 600000
        }
    )

#write the input file and run the simulation
from pyalps.run_io import execute, write_run_files
runs = [dict(parameters=p, output=dict(results=f'parm5b.task{i+1}.out.h5'))
        for i, p in enumerate(parms)]
files = execute('worm', write_run_files('parm5b', runs, baseseed=42, overwrite=True))

#load the magnetization and collect it as function of field h
data = pyalps.loadMeasurements(files,'Stiffness')
rhos = pyalps.collectXY(data,x='t',y='Stiffness',foreach=['L'])

# multiply with the system size for the scaling plot
for s in rhos:
  s.y = s.y * float(s.props['L'])

#make plot
plt.figure()
pyalps.plot.plot(rhos)
plt.xlabel('Hopping $t/U$')
plt.ylabel('$\\rho _sL$')
plt.legend()
plt.title('Scaling plot for Bose-Hubbard model')
plt.show()
