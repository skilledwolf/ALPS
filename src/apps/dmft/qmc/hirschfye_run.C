// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#include "hirschfyesim.h"
#include "run_config.h"
#include "parallel.hpp"
#include <alps/alea/hdf5.hpp>
#include <alps/ngs/signal.hpp>
#include <chrono>
#include <optional>
#include <stdexcept>
#include <string>
#ifdef ALPS_HAVE_MPI
#include <alps/check_schedule.hpp>
#endif

// Every rank calls at the same stage/check. A failed chain waits here while
// healthy chains reach their next scheduled check; no per-step collective.
void alps::dmft::agree_hirschfye_failure(std::string const& failure) {
  alps::solvers::parallel{}.agree(failure);
}

void alps::dmft::run_hirschfye(run_configuration run) {
  std::string failure;
  alps::solvers::parallel group;
  const int rank=group.rank;
  std::optional<HirschFyeRun> owned;
  try {
    const auto g0=prepare_hirschfye_run(run);
    owned.emplace(run,g0,rank);
  } catch (std::exception const& error) { failure=error.what(); }
  agree_hirschfye_failure(failure);
  auto& simulation=*owned;
  const auto started=std::chrono::steady_clock::now();
  const auto seconds=run.execution["time_limit"].as<int>();
  alps::ngs::signal signal;
  const auto stop=[&] {
    return !signal.empty() || (seconds && std::chrono::steady_clock::now()-started >=
                                        std::chrono::seconds(seconds));
  };
#ifdef ALPS_HAVE_MPI
  if (group.size==1) { while (!stop() && simulation.work_done()<1.) simulation.dostep(); }
  else {
  alps::check_schedule check;
  double fraction=0.;
  do {
    try { simulation.dostep(); }
    catch (std::exception const& error) { failure=error.what(); }
    if (!failure.empty() || check.pending()) {
      agree_hirschfye_failure(failure);
      fraction=stop() ? 1. : simulation.work_done();
      MPI_Allreduce(MPI_IN_PLACE,&fraction,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      check.update(fraction);
    }
  } while (fraction<1.);
  }
#else
  while (!stop() && simulation.work_done()<1.) simulation.dostep();
#endif
  try {
    const auto results=group.collect(simulation);
    if (!rank) {
      hdf5::save_checkpoint(run.output["results"].as<std::string>(),[&](hdf5::archive& archive) {
        archive["/parameters"]<<run.parameters;
        archive.create_group("/simulation/results");
        alea::hdf5_serializer serializer(archive,"/simulation/results");
        for (auto const& [name,result] : results)
          if (result.count()) serialize(serializer,archive.encode_segment(name),result);
        if (results.at("Sign").count()) {
          const auto green=HirschFyeRun::get_result(results,run.parameters);
          green.first.write_hdf5(archive,"/G_omega");
          green.second.write_hdf5(archive,"/G_tau");
        }
        archive["/run_config"]<<run;
      });
    }
  } catch (std::exception const& error) { failure=error.what(); }
  agree_hirschfye_failure(failure);
}
