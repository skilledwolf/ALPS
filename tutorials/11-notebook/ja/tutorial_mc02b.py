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
from pyalps.run_io import execute, write_run_files
import matplotlib.pyplot as plt
import pyalps.plot

#prepare the input parameters
parms = []
for t in [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0]:
    parms.append(
        { 
          'LATTICE'        : "ladder", 
          'T'              : t,
          'J0'             : [-1] ,
          'J1'             : [-1] ,
          'THERMALIZATION' : 10000,
          'SWEEPS'         : 500000,
          'UPDATE'         : "cluster",
          'MODEL'          : "Heisenberg",
          'L'              : 60
        }
    )

# Write typed task files and execute the job manifest.
input_file = write_run_files('mc02b', [{"parameters": p,
    "output": {"results": f"mc02b.task{i + 1}.out.h5"}}
    for i, p in enumerate(parms)], baseseed=42, overwrite=True)
execute('spinmc', input_file)

#load the susceptibility and collect it as function of temperature T
data = pyalps.loadMeasurements(pyalps.getResultFiles(prefix='mc02b'),'Susceptibility')
susceptibility = pyalps.collectXY(data,x='T',y='Susceptibility')

#make plot
plt.figure()
pyalps.plot.plot(susceptibility)
plt.xlabel('Temperature $T/J$')
plt.ylabel(r'Susceptibility $\chi J$')
plt.ylim(0,0.22)
plt.title('Heisenberg ladder')
plt.show()
