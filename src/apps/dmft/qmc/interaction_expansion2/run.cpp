/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2010 by Emanuel Gull <gull@phys.columbia.edu>,
 *
 *
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
 *
 *****************************************************************************/

#include <alps/solvers.hpp>
#include <alps/ctint.hpp>
#include "interaction_expansion.hpp"
#include <alps/utility/copyright.hpp>
#include <chrono>
#ifdef ALPS_HAVE_MPI
#include <alps/mcmpiadapter.hpp>
using sim_type = alps::mcmpiadapter<HubbardInteractionExpansion>;
#else
using sim_type = HubbardInteractionExpansion;
#endif

void compute_greens_functions(const alps::results_type<HubbardInteractionExpansion>::type&,
                             alps::params const&, alps::params const&, alps::params const&);

void alps::solvers::ctint(const run_configuration &supplied) {
  auto run = supplied;
  alps::ctint::prepare_run(run);
  int rank;
#ifndef ALPS_HAVE_MPI
  rank=0;
  sim_type s(run,rank);
#else
  boost::mpi::communicator c;
  c.barrier();
  rank=c.rank();
  const auto interval = run.execution["check_interval"].as<double>();
  sim_type s(run, c, alps::check_schedule(interval, interval),
             run.execution["bins"].as<std::size_t>());
#endif
  if (rank==0) {
    alps::print_copyright(std::cout);
    std::cout << "****************************************************************"<<std::endl;
    std::cout << "* Recommended citation in scientific publications:             *"<<std::endl;
    std::cout << "* We used the ALPS [1] implementation [2] of the CT-INT        *"<<std::endl;
    std::cout << "* interaction expansion CT-QMC [3,4] solver.                   *"<<std::endl;
    std::cout << "* [1] JSTAT (2011) P05001; [2] CPC 182, 1078 (2011);           *"<<std::endl;
    std::cout << "* [3] PRB 72, 035122 (2005); [4] RMP 83, 349 (2011).           *"<<std::endl;
    std::cout << "****************************************************************"<<std::endl;
  }
  const auto started = std::chrono::steady_clock::now();
  const auto seconds = run.execution["time_limit"].as<int>();
  alps::ngs::signal signal;
  s.run([&] {
    return !signal.empty() || (seconds != 0 && std::chrono::steady_clock::now() - started >=
                               std::chrono::seconds(seconds));
  });

  // All MPI ranks participate in collection; only root writes results.
  auto results = collect_results(s);
  if (rank==0) {
    const auto output_file = run.output["results"].as<std::string>();
    save_results(results, run.parameters, output_file, "/simulation/results");
    if (results["Sign"].count() != 0)
      compute_greens_functions(results, run.parameters, run.input, run.output);
    alps::hdf5::archive archive(output_file, "a");
    archive["/run_config"] << run;
  }
#ifdef ALPS_HAVE_MPI
  c.barrier();
#endif
}
