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
#include <alps/alea/hdf5.hpp>
#include "hyb.hpp"
#include "hybevaluate.hpp"
#include <alps/utility/copyright.hpp>
#include <chrono>
#include <optional>
#include "../parallel.hpp"
#ifdef ALPS_HAVE_MPI
#include <alps/check_schedule.hpp>
#endif


int global_mpi_rank;

namespace {
void master_final_tasks(hybridization::results_type const& results,
                        alps::run_configuration const& run, alps::hdf5::archive& solver_output){
  //do some post processing: collect Green functions and write
  //them into hdf5 files; calls compute vertex at the very end

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
  parallel group;
  std::optional<hybridization> owned;
  global_mpi_rank=group.rank;
  group.checked([&] { alps::cthyb::prepare_run(run); owned.emplace(run,group.rank); });
  const auto &parms=run.parameters;
  const auto output_file=run.output["results"].as<std::string>();
  auto& s=*owned;
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
  const auto started=std::chrono::steady_clock::now();
  const auto seconds=run.execution["time_limit"].as<int>();
  alps::ngs::signal signal;
  const auto stop=[&] {
    return !signal.empty() || (seconds != 0 && std::chrono::steady_clock::now()-started >=
                               std::chrono::seconds(seconds));
  };
#ifdef ALPS_HAVE_MPI
  if (group.size==1) { s.run(stop); }
  else {
  // SWEEPS is aggregate work across independent chains, as in the old driver.
  alps::check_schedule check;
  double fraction=0.;
  std::string failure;
  do {
    try { s.update(); s.measure(); }
    catch (std::exception const& error) { failure=error.what(); }
    if (!failure.empty() || check.pending()) {
      group.agree(failure);
      fraction=stop() ? 1. : s.fraction_completed();
      MPI_Allreduce(MPI_IN_PLACE, &fraction, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      check.update(fraction);
    }
  } while (fraction < 1.);
  }
#else
  s.run(stop);
#endif

  // Every rank participates in collection and receives any output error.
  std::string output_error;
  try {
    auto results=group.collect(s);
    if(global_mpi_rank==0){
      if(!results.at("Sign").count())
        throw std::runtime_error("CT-HYB stopped before any measurements; no results were written");
      const auto output_path=run.output["base_path"].as<std::string>()+"/simulation/results";
      alps::hdf5::save_checkpoint(output_file, [&](alps::hdf5::archive& output) {
        output["/parameters"]<<parms;
        output.create_group(output_path);
        alps::alea::hdf5_serializer serializer(output, output_path);
        for (auto const& entry : results)
          if (entry.second.count())
            serialize(serializer, output.encode_segment(entry.first), entry.second);
        master_final_tasks(results,run,output);
        output["/run_config"]<<run;
      });
    }
  } catch(const std::exception& error) { output_error=error.what(); }
  group.agree(output_error);
}
