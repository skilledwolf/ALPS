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
#include <alps/ngs/signal.hpp>
#include <alps/alea/hdf5.hpp>
#include "interaction_expansion.hpp"
#include "run_config.hpp"
#include <alps/utility/copyright.hpp>
#include <chrono>
#include <optional>
#ifdef ALPS_HAVE_MPI
#include <alps/alea/mpi.hpp>
#include <alps/check_schedule.hpp>
#include <climits>
#endif

void compute_greens_functions(InteractionExpansion::results_type const&,
                             alps::params const&, alps::params const&,
                             alps::hdf5::archive&);

void alps::ctint::agree_failure(std::string const& failure) {
#ifdef ALPS_HAVE_MPI
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  int origin = failure.empty() ? size : rank;
  MPI_Allreduce(MPI_IN_PLACE, &origin, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  if (origin == size) return;
  auto message = rank == origin ? failure : std::string{};
  if (message.size() > INT_MAX) message = "CT-INT error exceeds MPI message limit";
  int length = static_cast<int>(message.size());
  MPI_Bcast(&length, 1, MPI_INT, origin, MPI_COMM_WORLD);
  message.resize(length);
  MPI_Bcast(message.data(), length, MPI_CHAR, origin, MPI_COMM_WORLD);
  throw std::runtime_error(message);
#else
  if (!failure.empty()) throw std::runtime_error(failure);
#endif
}

void alps::solvers::ctint(const run_configuration &supplied) {
  auto run = supplied;
  int rank=0;
#ifdef ALPS_HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  std::string failure;
  std::optional<InteractionExpansion> owned;
  try {
    alps::ctint::prepare_run(run);
    owned.emplace(run, rank);
  } catch (std::exception const& error) { failure = error.what(); }
  alps::ctint::agree_failure(failure);
  auto& s = *owned;
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
  const auto stop = [&] {
    return !signal.empty() || (seconds != 0 && std::chrono::steady_clock::now() - started >=
                               std::chrono::seconds(seconds));
  };

#ifdef ALPS_HAVE_MPI
  int processes=1;
  MPI_Comm_size(MPI_COMM_WORLD, &processes);
  if (processes==1) { s.run(stop); }
  else {
  // SWEEPS is aggregate work across independent chains, as in the old driver.
  const auto interval = run.execution["check_interval"].as<double>();
  alps::check_schedule check(interval, interval);
  double fraction=0.;
  do {
    try { s.update(); s.measure(); }
    catch (std::exception const& error) { failure = error.what(); }
    if (!failure.empty() || check.pending()) {
      alps::ctint::agree_failure(failure);
      fraction = stop() ? 1. : s.fraction_completed();
      MPI_Allreduce(MPI_IN_PLACE, &fraction, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      check.update(fraction);
    }
  } while (fraction < 1.);
  }
  alps::alea::mpi_reducer reducer(MPI_COMM_WORLD);
#else
  s.run(stop);
#endif
  // All MPI ranks participate in collection; only root writes results.
  try {
#ifdef ALPS_HAVE_MPI
    auto results = s.collect_results(&reducer);
#else
    auto results = s.collect_results();
#endif
    if (rank==0) {
      const auto output_file = run.output["results"].as<std::string>();
      alps::hdf5::save_checkpoint(output_file, [&](alps::hdf5::archive& archive) {
        archive["/parameters"] << run.parameters;
        const std::string path = "/simulation/results";
        archive.create_group(path);
        alps::alea::hdf5_serializer serializer(archive, path);
        for (auto const& entry : results)
          if (entry.second.count())
            serialize(serializer, archive.encode_segment(entry.first), entry.second);
        archive["/run_config"] << run;
        if (results.at("Sign").count() != 0)
          compute_greens_functions(results, run.parameters, run.input, archive);
      });
    }
  } catch (const std::exception& error) { failure = error.what(); }
  alps::ctint::agree_failure(failure);
}
