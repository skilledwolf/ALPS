/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2002-2006 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

// How the legacy scheduler pools the measurements of a task's runs into the
// task's results. The merge was moved unchanged from MCSimulation
// (montecarlo.C): accumulate_measurements and the end of get_measurements.
//
// Runs are added in order and merged like a binary counter: whenever the two
// last groups hold the same number of runs, the later merges into the
// earlier. When the results are taken, the remaining groups merge from the
// right. Merging combines each observable into one and keeps the runs' data
// in run order; runs are not kept apart. For the scalar observables the tests
// cover, the grouping does not change the results but the run order does, so
// a backend that pools replicas in run order reproduces the legacy results.
//
// MCSimulation::get_and_remove_observable does not pool this way.

#ifndef ALPS_SCHEDULER_POOLING_H
#define ALPS_SCHEDULER_POOLING_H

#include <alps/alea/observableset.h>
#include <cstddef>
#include <utility>
#include <vector>

namespace alps {
namespace scheduler {

class measurement_pool
{
public:
  /// Adds the measurements of the next run.
  void add(const ObservableSet& measurements)
  {
    groups_.push_back(std::make_pair(std::size_t(1), measurements));
    while (groups_.size() > 1 && groups_.back().first == (groups_.rbegin() + 1)->first) {
      (groups_.rbegin() + 1)->first *= 2;
      (groups_.rbegin() + 1)->second << groups_.back().second;
      groups_.pop_back();
    }
  }

  bool empty() const { return groups_.empty(); }

  /// The pooled measurements of the runs added so far, with the observables
  /// of simulation whose names no run has, and empties the pool. An empty
  /// pool gives an empty set.
  ObservableSet take(const ObservableSet& simulation)
  {
    if (groups_.empty())
      return ObservableSet();
    for (std::size_t i = groups_.size() - 1; i > 0; --i)
      groups_[i - 1].second << groups_[i].second;
    for (ObservableSet::const_iterator it=simulation.begin();it!=simulation.end();++it)
      if (!groups_[0].second.has(it->first))
        groups_[0].second << *(it->second);
    ObservableSet pooled(groups_[0].second);
    groups_.clear();
    return pooled;
  }

private:
  std::vector<std::pair<std::size_t, ObservableSet> > groups_;
};

} // namespace scheduler
} // namespace alps

#endif // ALPS_SCHEDULER_POOLING_H
