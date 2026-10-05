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
for t in [1.5,2,2.5]:
    parms.append(
        { 
          'LATTICE'        : "square lattice", 
          'T'              : t,
          'J'              : [1] ,
          'THERMALIZATION' : 1000,
          'SWEEPS'         : 100000,
          'UPDATE'         : "cluster",
          'MODEL'          : "Ising",
          'L'              : 8
        }
    )

# Write typed task files and execute the job manifest.
input_file = write_run_files('parm1', [{"parameters": p,
    "output": {"results": f"parm1.task{i + 1}.out.h5"}}
    for i, p in enumerate(parms)], baseseed=42, overwrite=True)
execute('spinmc', input_file)

#get the list of result files
result_files = pyalps.getResultFiles(prefix='parm1')
print("Loading results from the files: ", result_files)

#print the observables stored in those files:
print("The files contain the following mesurements:", end=' ')
print(pyalps.loadObservableList(result_files))

# The native joint jackknife retains covariance between the moments.
data = pyalps.loadMeasurements(result_files, 'Binder Cumulant U2')
binder = pyalps.collectXY(data, 'T', 'Binder Cumulant U2')

# ... and plot them
plt.figure()
pyalps.plot.plot(binder)
plt.xlabel('T')
plt.ylabel('Binder cumulant')
plt.show()
