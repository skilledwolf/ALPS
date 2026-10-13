/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2010 by Matthias Troyer <troyer@comp-phys.org>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

// How the legacy scheduler measures the progress of a task, in one place. The
// expressions were moved unchanged from Task::work (task.C) and from
// WorkerTask::work and WorkerTask::finished (workertask.C).
//
// Each run reports the fraction of its work it has done. A task's work done
// is the sum over its runs, including runs on dump, so a budget such as
// SWEEPS counts for the whole task rather than for each run. The task is
// finished once its remaining work is no longer positive.

#ifndef ALPS_SCHEDULER_PROGRESS_H
#define ALPS_SCHEDULER_PROGRESS_H

#include <alps/expression.h>
#include <alps/parameter.h>

namespace alps {
namespace scheduler {

/// WORK_FACTOR evaluated as an expression, or 1 if it is not defined. The
/// work of a task without runs.
inline double work_factor(const Parameters& parms)
{
  return parms.defined("WORK_FACTOR") ? alps::evaluate<double>(parms["WORK_FACTOR"], parms) : 1.;
}

/// The work a task with runs still has to do: work_factor times the fraction
/// of work not yet done, which turns negative once the runs overshoot.
inline double remaining_work(const Parameters& parms, double work_done)
{
  return work_factor(parms)*(1.-work_done);
}

/// The percentage reported for a task with the given remaining work: 1 minus
/// that work, clipped to [0, 1]. With a WORK_FACTOR other than 1 it is not
/// the fraction of work done.
inline double reported_percentage(double work)
{
  double percentage = 1.-work;
  if (percentage < 0.)
    percentage=0.;
  else if (percentage>1.)
    percentage=1.;
  return percentage;
}

} // namespace scheduler
} // namespace alps

#endif // ALPS_SCHEDULER_PROGRESS_H
