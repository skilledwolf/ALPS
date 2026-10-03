/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2009 by Emanuel Gull <gull@phys.columbia.edu>
 *                              Philipp Werner <werner@itp.phys.ethz.ch>,
 *                              Sebastian Fuchs <fuchs@theorie.physik.uni-goettingen.de>
 *                              Matthias Troyer <troyer@comp-phys.org>
 *
 *
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
 *
 *****************************************************************************/

/* $Id: alps_solver.C 342 2009-01-28 22:31:54Z fuchs $ */

#include "alps_solver.h"
#include "xml.h"
#include "types.h"
#include <boost/assert.hpp>
#include <alps/osiris/comm.h>
#include <alps/scheduler/options.h>
#include <alps/ngs/make_deprecated_parameters.hpp>
#include <limits>
#include <vector>
#include <utility>
#include <sstream>

#define DEFAULT_CHECK_TIME 300

alps::ImpuritySolver::ImpuritySolver(const scheduler::Factory& factory, const run_configuration& run, int argc, char** argv)
  : master_scheduler(nullptr), configuration_(run)
{
#ifdef ALPS_HAVE_MPI
  comm_init(argc,argv, true);
#else
  comm_init(argc,argv, false);
#endif
  alps::scheduler::NoJobfileOptions opt(1,argv);
  const auto seconds = run.execution["time_limit"].as<unsigned int>();
  const auto interval = seconds ? std::min(seconds, unsigned(DEFAULT_CHECK_TIME)) : unsigned(DEFAULT_CHECK_TIME);
  opt.checkpoint_time = opt.max_check_time = opt.min_check_time = interval;
  if (is_master()) master_scheduler = new alps::scheduler::SingleScheduler(opt,factory);
  else {
    alps::scheduler::Scheduler slave(opt,factory);
    slave.run(); comm_exit(); exit(0);
  }
}


alps::ImpuritySolver::~ImpuritySolver()
{
  scheduler::stop_single();
}


int alps::ImpuritySolver::solve_it(Parameters const& p)
{
  if (is_master())
  {
    master_scheduler->create_task(p);
    return master_scheduler->run();
  }
  else{
    std::cerr<<"a slave should never reach this section."<<std::endl;
    abort();
  }
}


std::pair<matsubara_green_function_t, itime_green_function_t>
alps::ImpuritySolver::solve_omega(const matsubara_green_function_t& G0_omega, const params& p)
{
  BOOST_ASSERT(is_master());
  // The scheduler tasks still read untyped Parameters; U_MATRIX names the
  // interaction file again because MCRun cannot receive the typed input section.
  alps::Parameters parms = make_deprecated_parameters(p);
  parms["SEED"] = configuration_.execution["seed"].as<std::uint64_t>();
  const auto seconds = configuration_.execution["time_limit"].as<int>();
  parms["MAX_TIME"] = seconds ? seconds : std::numeric_limits<int>::max();
  parms["OUTFILE"] = configuration_.output["results"].as<std::string>();
  if (configuration_.input.exists("interaction_matrix"))
    parms["U_MATRIX"] = configuration_.input["interaction_matrix"].as<std::string>();
  
  std::ostringstream G0_omega_text;
  alps::oxstream G0_omega_xml(G0_omega_text);
  write_freq(G0_omega_xml,G0_omega);
  parms["G0(omega)"] = G0_omega_text.str();
  if (solve_it(parms))
    boost::throw_exception(std::runtime_error(" solver finished with nonzero exit code"));
  
  std::pair<matsubara_green_function_t, itime_green_function_t> G = 
    dynamic_cast<MatsubaraImpurityTask*>(get_task())->get_result();
  clear(); // destroy the simulation
  return G;
}
