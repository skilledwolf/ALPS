/****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2012 by Emanuel Gull <gull@pks.mpg.de>,
 *                       Hartmut Hafermann <hafermann@cpht.polytechnique.fr>
 *
 *  based on an earlier version by Philipp Werner and Emanuel Gull
 *
 *
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
 *
 *****************************************************************************/

#include <alps/solvers.hpp>
#include <alps/cthyb.hpp>
#include <alps/ngs/signal.hpp>
#include "hyb.hpp"
#include "hybevaluate.hpp"
#include <alps/utility/copyright.hpp>
#include <boost/date_time/posix_time/posix_time_types.hpp>
#ifdef ALPS_HAVE_MPI
#include <alps/mcmpiadapter.hpp>
using sim_type = alps::mcmpiadapter<hybridization>;
#else
using sim_type = hybridization;
#endif


int global_mpi_rank;

namespace {
bool stop_requested(boost::posix_time::ptime const & end_time) {
  static alps::ngs::signal signal;
  return !signal.empty() || (!end_time.is_not_a_date_time() && boost::posix_time::second_clock::local_time() > end_time);
}
void master_final_tasks(const alps::results_type<hybridization>::type &results,
                        const alps::run_configuration &run){
  //do some post processing: collect Green functions and write
  //them into hdf5 files; calls compute vertex at the very end

  alps::hdf5::archive solver_output(run.output["results"].as<std::string>(), "a");
  const auto &parms=run.parameters;

  evaluate_basics(results,parms,run.output,solver_output);
  evaluate_time(results,parms,run.output,solver_output);
  evaluate_freq(results,parms,run.output,solver_output);
  evaluate_legendre(results,parms,run.output,solver_output);
  evaluate_nnt(results,parms,run.output,solver_output);
  evaluate_nnw(results,parms,run.output,solver_output);
  evaluate_sector_statistics(results,parms,run.output,solver_output);
  evaluate_2p(results, parms, run.output, solver_output);
}

}

void alps::solvers::cthyb(alps::run_configuration const& supplied) {
  auto run=supplied;
  alps::cthyb::prepare_run(run);
  const auto &parms=run.parameters;
  const auto output_file=run.output["results"].as<std::string>();
#ifndef ALPS_HAVE_MPI
  global_mpi_rank=0;
  sim_type s(run,global_mpi_rank);
#else
  boost::mpi::communicator c;
  c.barrier();
  global_mpi_rank=c.rank();
  sim_type s(run, c, alps::check_schedule(), run.execution["bins"].as<std::size_t>());
#endif
  if (global_mpi_rank==0) {
    alps::print_copyright(std::cout);
    std::cout << "****************************************************************"<<std::endl;
    std::cout << "* Recommended citation in scientific publications:             *"<<std::endl;
    std::cout << "* We used the ALPS [1] implementation [2] of the `segment'     *"<<std::endl;
    std::cout << "* CT-QMC solver [3,4,5].                                       *"<<std::endl;
    std::cout << "* [1] JSTAT (2011) P05001; [2] CPC 182, 1078 (2011);           *"<<std::endl;
    std::cout << "* [3] PRL 97, 076405 (2006); [4] RMP 83, 349 (2011);           *"<<std::endl;
    std::cout << "* [5] PRB 84, 075145 (2011).                                   *"<<std::endl;
    std::cout << "****************************************************************"<<std::endl;
  }
  //run the simulation
  const auto time_limit=run.execution["time_limit"].as<int>();
  const auto end=time_limit ? boost::posix_time::second_clock::local_time()+boost::posix_time::seconds(time_limit)
                            : boost::posix_time::ptime(boost::posix_time::not_a_date_time);
  s.run(boost::bind(&stop_requested,end));

  // Every rank participates in collection and receives any output error.
  auto results=collect_results(s);
  std::string output_error;
  if(global_mpi_rank==0){
    try {
      if(!results["Sign"].count())
        throw std::runtime_error("CT-HYB stopped before any measurements; no results were written");
      const auto output_path=run.output["base_path"].as<std::string>()+"/simulation/results";
      save_results(results,parms,output_file,output_path);
      master_final_tasks(results,run);
      alps::hdf5::archive output(output_file,"a");
      output["/run_config"]<<run;
    } catch(const std::exception& error) { output_error=error.what(); }
  }
#ifdef ALPS_HAVE_MPI
  boost::mpi::broadcast(c,output_error,0);
#endif
  if(!output_error.empty()) throw std::runtime_error(output_error);
}
