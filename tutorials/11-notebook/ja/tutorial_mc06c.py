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
parms = [{ 
          'LATTICE'        : "simple cubic lattice", 
          'MODEL'          : "spin",
          'local_S'        : 0.5,
          'L'              : 4,
          'J'              : 1 ,
          'CUTOFF'         : 500
        }]

#write the input file and run the simulation
runs = [{"parameters": p, "execution": {"seed": 42 + i},
         "output": {"results": f"mc06c.task{i + 1}.out.h5"}} for i, p in enumerate(parms)]
files = execute("qwl", write_run_files("mc06c", runs, overwrite=True))

#run the evaluation and load all the plots
data = pyalps.evaluateQWL(files,DELTA_T=0.05, T_MIN=0.5, T_MAX=5.0)

#make plot
for s in pyalps.flatten(data):
  plt.figure()
  plt.title("Cubic lattice Heisenberg antiferromagnet L=4")
  pyalps.plot.plot(s)

plt.show()
