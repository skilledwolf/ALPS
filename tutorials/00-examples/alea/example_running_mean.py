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

# This is an example of how to calculate the running mean and the reverse running mean of data stored in a hdf5 file.

# Run generate_samples.py first. These are raw observations, not bin means.
filename = "timeseries.h5"
with pyalps.hdf5.archive(filename) as archive:
    obs = archive["/samples/E"]

# calculate the running mean
running_mean = pyalps.alea.running_mean(obs)

# calculate the reverse running mean
reverse_running_mean = pyalps.alea.reverse_running_mean(obs)

# print the result
print("The running mean of E is: " + str(running_mean))
print()
print("The reverse running mean of E is: " + str(reverse_running_mean))
print()


