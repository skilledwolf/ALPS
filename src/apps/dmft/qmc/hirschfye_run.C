// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#include "hirschfyesim.h"
#include "run_config.h"
#include <alps/alea/hdf5.hpp>
#include <alps/ngs/signal.hpp>
#include <chrono>
#include <optional>
#include <stdexcept>
#include <string>
#ifdef ALPS_HAVE_MPI
#include <alps/alea/mpi.hpp>
#include <alps/check_schedule.hpp>
#include <climits>
#endif

// Every rank calls at the same stage/check. A failed chain waits here while
// healthy chains reach their next scheduled check; no per-step collective.
void alps::dmft::agree_hirschfye_failure(std::string const& failure) {
#ifdef ALPS_HAVE_MPI
  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  int origin=failure.empty() ? size : rank;
  MPI_Allreduce(MPI_IN_PLACE,&origin,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  if (origin==size) return;
  auto message=rank==origin ? failure : std::string{};
  if (message.size()>INT_MAX) message="Hirsch-Fye error exceeds MPI message limit";
  int length=static_cast<int>(message.size());
  MPI_Bcast(&length,1,MPI_INT,origin,MPI_COMM_WORLD);
  message.resize(length);
  MPI_Bcast(message.data(),length,MPI_CHAR,origin,MPI_COMM_WORLD);
  throw std::runtime_error(message);
#else
  if (!failure.empty()) throw std::runtime_error(failure);
#endif
}

void alps::dmft::run_hirschfye(run_configuration run) {
  std::string failure;
  int rank=0;
#ifdef ALPS_HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
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
  alps::alea::mpi_reducer reducer(MPI_COMM_WORLD);
#else
  while (!stop() && simulation.work_done()<1.) simulation.dostep();
#endif
  try {
#ifdef ALPS_HAVE_MPI
    const auto results=simulation.collect_results(&reducer);
#else
    const auto results=simulation.collect_results();
#endif
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
