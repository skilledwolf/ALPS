# ****************************************************************************
#
# ALPS Project: Algorithms and Libraries for Physics Simulations
#
# ALPS Libraries
#
# Copyright (C) 2010 by Bela Bauer <bauerb@phys.ethz.ch>
#
# ALPS Project: https://alps.comp-phys.org/
# SPDX-License-Identifier: MIT
#
# ****************************************************************************
"""Seeds that writeInputFiles gives to tasks.

The derivation was moved unchanged out of writeInputFiles (tools.py). It has
a C++ counterpart, alps::scheduler::task_seed in alps/scheduler/seeds.h, which
also describes how the scheduler derives the seeds of runs from a task's SEED.
Both are checked against src/alps/scheduler/tests/seed_vectors.txt.

This module imports nothing from pyalps, so it can be used and tested without
the compiled extensions.
"""


def task_seed(baseseed, index, ntasks):
    """Return the SEED for task ``index`` (counted from 0) of ``ntasks``.

    The index is folded into ``baseseed`` with XOR, repeated every ``bits``
    bits, where ``bits`` shrinks as the number of tasks grows. The result is
    masked to 31 bits.
    """
    bits = 31
    n = ntasks
    while n>0:
      n //= 2
      bits -= 1

    seed = baseseed
    for j in range(0,32//bits+1):
      seed ^= (index << (j * bits))
      seed &= ((1<<30) | ((1<<30)-1))
    return seed
