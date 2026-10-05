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
import pyalps.plot
import subprocess

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

#load a selection of measurements:
data = pyalps.loadMeasurements(result_files,['|Magnetization|','Magnetization^2'])

#make a plot for the magnetization: collect Magnetziation as function of T
plotdata = pyalps.collectXY(data,'T','|Magnetization|')

# convert the data to grace file for plotting using xmgrace
print("The results in grace format are:")
print(pyalps.plot.makeGracePlot(plotdata))

print("Saving into file parm1.gr")
f = open ('parm1.gr','w')
f.write(pyalps.plot.makeGracePlot(plotdata))
f.close()

print("Trying to launch xmgrace")
process = subprocess.Popen(['xmgrace','parm1.gr'])

