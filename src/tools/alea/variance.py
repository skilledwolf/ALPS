#! /usr/bin/env python
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

from mcanalyze_tools import impl_calculation

def calculate (obs):
  if not hasattr(obs, "variance"):
    raise ValueError("This ALEA result stores a mean only; variance cannot be recovered")
  return obs.variance

impl_calculation("Variance", "variance/value", calculate)

