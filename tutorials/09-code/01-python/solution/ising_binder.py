# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2009-2010 by Matthias Troyer <troyer@phys.ethz.ch> 
#                            Jan Gukelberger
#                            Brigitte Surer
# 
# ALPS Project: https://alps.comp-phys.org/
# SPDX-License-Identifier: MIT
# 
# ****************************************************************************

# The same simulation retains m^2/m^4 covariance in its joint ALEA batches.
from ising import Simulation, main

if __name__ == '__main__':
    main(sizes=(4, 6, 8), plot=True)
