/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2010 by Matthias Troyer <troyer@comp-phys.org>,
*                            Synge Todo <wistaria@comp-phys.org>
* Copyright (C) 2010 by Bela Bauer <bauerb@phys.ethz.ch>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

// How legacy ALPS derives random number seeds, in one place. The expressions
// were moved unchanged from the scheduler: worker_seed and disorder_seed from
// Worker (worker.C), next_run_seed from WorkerTask (workertask.C). task_seed
// restates the derivation in pyalps.writeInputFiles (tools.py, now
// pyalps/seeds.py) so a C++ front end can reproduce the seeds pyalps writes.
// src/alps/scheduler/tests/seed_vectors.txt holds reference values that the
// tests of both languages check.
//
// The seeds of one task, from the outside in:
//   1. pyalps gives task i of n the seed task_seed(baseseed, i, n), unless
//      the task already sets SEED. The task scheduler sets SEED=0 if absent.
//   2. The task gives its first run SEED and every further run the previous
//      value advanced by next_run_seed.
//   3. Each run seeds its engine with worker_seed and the process-wide
//      disorder generator with disorder_seed (DISORDERSEED, default 0, the
//      same for all runs).
// Parameter expressions also reseed the disorder generator whenever
// DISORDERSEED, read as disorder_seed reads it, is not 0 and differs from
// its last seed (alps::Disorder::seed_if_unseeded).
// Parapack and the ngs scheduler derive their seeds separately.

#ifndef ALPS_SCHEDULER_SEEDS_H
#define ALPS_SCHEDULER_SEEDS_H

#include <alps/parameter.h>
#include <boost/cstdint.hpp>
#include <boost/throw_exception.hpp>
#include <stdexcept>

namespace alps {
namespace scheduler {

/// Seed of a worker's random number engine: SEED read as uint32_t, so a
/// negative value wraps around. A missing SEED fails like an empty one,
/// with boost::bad_lexical_cast.
inline uint32_t worker_seed(const Parameters& parms)
{
  return parms.value_or_default("SEED","");
}

/// Seed of the process-wide disorder generator: DISORDERSEED read as
/// unsigned int, or 0.
inline unsigned int disorder_seed(const Parameters& parms)
{
  return parms.value_or_default("DISORDERSEED",0);
}

/// SEED of the next run of a task, after a run with SEED seed that uses
/// cpus processes. The sum int32_t(seed) + cpus is unsigned, so a negative
/// SEED advances to a large value. A SEED above INT32_MAX cannot be read
/// and fails with boost::bad_lexical_cast.
inline uint32_t next_run_seed(const StringValue& seed, uint32_t cpus)
{
  return static_cast<int32_t>(seed)+cpus;
}

/// SEED that pyalps.writeInputFiles gives task index (from 0) of ntasks:
/// baseseed XOR the index repeated every bits bits, where bits shrinks with
/// the number of tasks, masked to 31 bits.
inline uint32_t task_seed(uint64_t baseseed, uint64_t index, uint64_t ntasks)
{
  int bits = 31;
  uint64_t n = ntasks;
  while (n>0) {
    n /= 2;
    bits -= 1;
  }
  if (bits<=0)
    boost::throw_exception(std::invalid_argument("too many tasks to derive seeds"));
  uint64_t seed = baseseed;
  for (int j=0; j<32/bits+1; ++j) {
    seed ^= (index << (j * bits));
    seed &= ((1u<<30) | ((1u<<30)-1));
  }
  return static_cast<uint32_t>(seed);
}

} // end namespace scheduler
} // end namespace alps

#endif // ALPS_SCHEDULER_SEEDS_H
