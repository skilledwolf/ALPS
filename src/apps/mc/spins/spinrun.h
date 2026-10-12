/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1999-2009 by Matthias Troyer <troyer@comp-phys.org>,
*                            Fabian Stoeckli <fabstoec@student.ethz.ch>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

// Runs the classical spin simulation under the legacy scheduler. Until the
// simulation was separated from the scheduler, AbstractSpinSim was itself a
// scheduler run (alps::scheduler::LatticeMCRun); this adapter keeps that
// run's behavior. It hands the simulation the worker's parameters,
// measurements and random numbers, forwards the scheduler's calls, and keeps
// what only concerns a scheduler run: a run without processes compacts its
// measurements and drops its spins when it is loaded (moved unchanged from
// AbstractSpinSim::load and SpinSim::load).

#ifndef ALPS_APPLICATIONS_MC_SPIN_SPINRUN_H_
#define ALPS_APPLICATIONS_MC_SPIN_SPINRUN_H_

#include <alps/scheduler.h>
#include "spinsim.h"

template <class M, class MAT>
class SpinRun : public alps::scheduler::MCRun {
public:
  SpinRun(const alps::ProcessList& w, const alps::Parameters& p, int n)
    : alps::scheduler::MCRun(w,p,n),
      sim_(parms,measurements,random_01)
  {}

  void dostep() override { sim_.dostep(); }
  double work_done() const override { return sim_.work_done(); }
  bool is_thermalized() const override { return sim_.is_thermalized(); }
  bool change_parameter(const std::string& name, const alps::StringValue& value) override
  { return sim_.change_parameter(name,value); }
  void save(alps::ODump& dump) const override { sim_.save(dump); }
  void load(alps::IDump& dump) override
  {
    sim_.load(dump);
    if (where.empty()) {
      measurements.compact();
      sim_.clear_spins();
    }
  }

private:
  SpinSim<M,MAT> sim_;
};

#endif
