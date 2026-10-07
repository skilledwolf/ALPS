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
parms = [{ 
          'LATTICE'                   : "double dimer", 
          'MODEL'                     : "dimerized spin",
          'CONSERVED_QUANTUMNUMBERS'  : 'Sz',
          'local_S0'                  : 1,
          'local_S1'                  : 0.5,
          'J0'                        : 1,
          'J1'                        : 0.4
        }]

#write the input file and run the simulation
from pyalps.run_io import execute, write_run_files
runs = [dict(parameters=p, input=dict(lattice_library='dd-graph.xml', model_library='model-dspin.xml'),
             output=dict(results=f'parm2c.task{i+1}.out.h5'))
        for i, p in enumerate(parms)]
files = execute('fulldiag', write_run_files('parm2c', runs, overwrite=True))

#run the evaluation and load all the plots
data = pyalps.evaluateFulldiagVersusH(files,T = 0.02, DELTA_H=0.025, H_MIN=0., H_MAX=4.0)

#make plot
for s in pyalps.flatten(data):
  plt.figure()
  pyalps.plot.plot(s)

plt.show()
