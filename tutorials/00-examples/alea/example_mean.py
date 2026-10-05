#/*****************************************************************************
#*
#* ALPS Project: Algorithms and Libraries for Physics Simulations
#*
#* Copyright (C) 2011-2012 by Lukas Gamper <gamperl@gmail.com>,
#*                            Matthias Troyer <troyer@itp.phys.ethz.ch>,
#*                            Maximilian Poprawe <poprawem@ethz.ch>
#*
#* ALPS Project: https://alps.comp-phys.org/
#* SPDX-License-Identifier: MIT
#*
#*****************************************************************************/

import pyalps

# This is an example of how to easily calculate the mean of data stored in a hdf5 file.


# Run generate_samples.py first. These are raw observations, not bin means.
filename = "timeseries.h5"
with pyalps.hdf5.archive(filename) as archive:
    obs = archive["/samples/E"]

# calculate the mean
mean = pyalps.alea.mean(obs)

# print the result
print("The mean of E is: " + str(mean))
  
# write the result back to the file
with pyalps.hdf5.archive(filename, "a") as archive:
    archive["/analysis/mean"] = mean
