// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#include "interaction_expansion.hpp"
#include <alps/ctint.hpp>
#include <alps/alea/mpi.hpp>
#include <alps/alea/hdf5.hpp>
#include <alps/hdf5/vector.hpp>
#include <alps/hdf5/complex.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>

namespace {
void require(bool value, char const* message) {
  int agreed = value;
  MPI_Allreduce(MPI_IN_PLACE, &agreed, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  if (!agreed) throw std::runtime_error(message);
}
template<class F> void rejects(F operation) {
  bool failed = false;
  try { operation(); } catch (std::exception const&) { failed = true; }
  require(failed, "invalid CT-INT collection was not rejected on every rank");
}
alps::run_configuration configuration() {
  alps::run_configuration run;
  run.parameters["BETA"] = 2.;
  run.parameters["U"] = 0.;
  run.parameters["MU"] = 0.;
  run.parameters["ALPHA"] = -.01;
  run.parameters["N"] = 8;
  run.parameters["NMATSUBARA"] = 4;
  run.parameters["SWEEPS"] = 8;
  run.parameters["THERMALIZATION"] = 2;
  run.parameters["MEASUREMENT_PERIOD"] = 1;
  run.input["atomic"] = true;
  run.output["results"] = "ctint-mpi-contract.h5";
  run.execution["bins"] = 16;
  run.execution["check_interval"] = 1.e-7;
  alps::ctint::prepare_run(run);
  return run;
}
struct fixture : HubbardInteractionExpansion {
  using HubbardInteractionExpansion::HubbardInteractionExpansion;
  void sample(double s, double value) {
    sign = s;
    record_measurement("Sign", s);
    record_measurement("densities", std::valarray<double>{s * value, 2 * s * value});
  }
  void change_name() {
    auto entry = measurements.extract("VertexRemoval");
    entry.key() = "VertexRemovas"; // Same registry count and name length.
    measurements.insert(std::move(entry));
  }
  void change_shape() { measurements.at("densities").accumulator.set_size(4); }
  void change_sign() { measurements.at("densities").signed_value = false; }
};
void collection(alps::run_configuration const& run, MPI_Comm comm) {
  int rank;
  MPI_Comm_rank(comm, &rank);
  alps::alea::mpi_reducer reducer(comm);
  for (int active : {0, 1, -1}) {
    fixture sim(run, rank);
    if (rank == active) for (int i=0; i<3; ++i) sim.sample(1., 2.);
    auto results = sim.collect_results(&reducer);
    bool ok = results.empty();
    if (!rank) {
      auto const& density = results.at("densities");
      ok = density.size() == 2 && density.count() == (active < 0 ? 0 : 3);
      if (active >= 0) ok = ok && density.mean()(0) == 2. && density.stderror()(0) == 0.;
    }
    require(ok, "empty CT-INT replica changed the root result");
  }
  fixture unequal(run, rank);
  for (int i=0; i<(rank ? 5 : 3); ++i) unequal.sample(rank && i>=3 ? -1. : 1., rank ? 5. : 2.);
  auto results = unequal.collect_results(&reducer);
  bool ok = results.empty();
  if (!rank) {
    auto const& density = results.at("densities");
    // Raw pooled ratio is 11/4. Unit-weight leave-one-out pseudovalues are
    // 1 (three), 8 (three), -.4 (two), giving mean 3.275 and error^2 1.955625.
    ok = results.at("Sign").count() == 8 && results.at("Sign").mean()(0) == .5
      && density.count() == 8 && std::abs(density.mean()(0)-3.275) < 1.e-12
      && std::abs(density.mean()(1)-6.55) < 1.e-12
      && std::abs(density.stderror()(0)-std::sqrt(1.955625)) < 1.e-12;
  }
  require(ok, "CT-INT reduced corrected replica means instead of joint signed batches");

  fixture zero_sign(run, rank);
  zero_sign.sample(rank ? -1. : 1., 2.);
  auto original = zero_sign.collect_results();
  rejects([&] { zero_sign.collect_results(&reducer); });
  require(zero_sign.collect_results() == original, "failed CT-INT analysis changed live measurements");
  for (int mismatch=0; mismatch<3; ++mismatch) {
    fixture invalid(run, rank);
    if (rank) {
      if (!mismatch) invalid.change_name();
      else if (mismatch == 1) invalid.change_shape();
      else invalid.change_sign();
    }
    rejects([&] { invalid.collect_results(&reducer); });
  }
}
void driver(alps::run_configuration run, int rank) {
  auto const path = run.output["results"].as<std::string>();
  require(rank || !std::filesystem::exists(path), "CT-INT test output already exists");
  alps::solvers::ctint(run);
  bool ok = true;
  std::string previous;
  if (!rank) try {
    alps::hdf5::archive archive(path);
    alps::alea::hdf5_serializer serializer(archive, "/simulation/results");
    alps::alea::batch_result<double> sign;
    deserialize(serializer, "Sign", sign);
    std::vector<std::complex<double>> green;
    archive["/G_omega/0/mean/value"] >> green;
    std::string application;
    archive["/run_config/application"] >> application;
    ok = sign.count() > 0 && sign.mean()(0) == 1. && application == "ctint"
      && green.size() == 4 && std::abs(green[0]-std::complex<double>(0., -2./std::acos(-1.))) < 1.e-12;
    archive.close();
    std::ifstream file(path, std::ios::binary);
    previous.assign(std::istreambuf_iterator<char>(file), {});
  } catch (std::exception const& error) { std::cerr << error.what() << '\n'; ok = false; }
  require(ok, "native CT-INT driver did not publish canonical results, Green functions and provenance");
  // Failure occurs inside Green-function processing, after result serialization.
  run.output["matrix_size"] = "ctint-mpi-missing-directory/matrix_size";
  rejects([&] { alps::solvers::ctint(run); });
  if (!rank) {
    std::ifstream file(path, std::ios::binary);
    ok = previous == std::string(std::istreambuf_iterator<char>(file), {});
    std::filesystem::remove(path);
  }
  require(ok, "failed CT-INT publication changed the previous output");
}
}
int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  int rank, size, status=0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  try {
    require(size == 2, "run CT-INT contract with two ranks");
    auto run = configuration();
    collection(run, MPI_COMM_WORLD);
    MPI_Comm reversed;
    MPI_Comm_split(MPI_COMM_WORLD, 0, 1-rank, &reversed);
    collection(run, reversed);
    MPI_Comm_free(&reversed);
    driver(run, rank);
  } catch (std::exception const& error) {
    if (!rank) std::cerr << error.what() << '\n';
    status=1;
  }
  MPI_Finalize();
  return status;
}
