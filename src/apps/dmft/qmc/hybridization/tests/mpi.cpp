// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#include "hyb.hpp"
#include <alps/cthyb.hpp>
#include <alps/solvers.hpp>
#include <alps/alea/mpi.hpp>
#include <alps/alea/hdf5.hpp>
#include <alps/hdf5/vector.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>

namespace {
void require(bool value, char const* message) {
  int agreed=value;
  MPI_Allreduce(MPI_IN_PLACE, &agreed, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  if (!agreed) throw std::runtime_error(message);
}
template<class F> void rejects(F operation) {
  bool failed=false;
  try { operation(); } catch (std::exception const&) { failed=true; }
  require(failed, "invalid CT-HYB operation was not rejected on every rank");
}
alps::run_configuration configuration() {
  alps::run_configuration run;
  run.parameters["BETA"]=2.;
  run.parameters["U"]=0.;
  run.parameters["N_ORBITALS"]=2;
  run.parameters["N_TAU"]=2;
  run.parameters["N_MEAS"]=2;
  // Complete at the initial MPI check, retaining the default check schedule.
  run.parameters["SWEEPS"]=1;
  run.parameters["THERMALIZATION"]=0;
  run.parameters["MEASURE_g2w"]=true;
  run.parameters["N_w2"]=2;
  run.parameters["N_W"]=1;
  run.input["delta"]="cthyb-mpi-delta.dat";
  run.output["results"]="cthyb-mpi-contract.h5";
  run.output["base_path"]="/impurity";
  run.execution["bins"]=16;
  alps::cthyb::prepare_run(run);
  return run;
}
struct fixture : hybridization {
  using hybridization::hybridization;
  using large_acc=alps::alea::var_acc<std::complex<double>,alps::alea::elliptic_var>;
  void sample(double inner_sign, double final_sign, double value) {
    record_measurement("Sign", inner_sign);
    record_measurement("density_0", inner_sign*value, inner_sign);
    record_measurement("g2w_re_0_0",
      std::vector<double>{final_sign*value, 2*final_sign*value,
                          3*final_sign*value, 4*final_sign*value}, final_sign);
  }
  void change_name() {
    auto entry=measurements.extract("order_0");
    entry.key()="order_9"; // Same registry count and name length.
    measurements.insert(std::move(entry));
  }
  void change_shape() {
    std::get<alps::alea::batch_acc<double>>(measurements.at("density_0").accumulator).set_size(3);
  }
  void change_family() { measurements.at("density_0").accumulator=large_acc(1); }
  void change_sign() { measurements.at("density_0").signed_value=false; }
  bool unchanged(double sign) const {
    auto raw=std::get<alps::alea::batch_acc<double>>(measurements.at("density_0").accumulator).result();
    auto large=std::get<large_acc>(measurements.at("g2w_re_0_0").accumulator).result();
    return raw.count()==1 && raw.mean()(0)==2*sign && raw.mean()(1)==sign
      && large.count()==1 && large.mean()(0)==std::complex<double>(2*sign,sign);
  }
};
void collection(alps::run_configuration const& run, MPI_Comm comm) {
  int rank;
  MPI_Comm_rank(comm, &rank);
  alps::alea::mpi_reducer reducer(comm);
  for (int active : {0, 1, -1}) {
    fixture sim(run,rank);
    if (rank==active) for (int i=0; i<3; ++i) sim.sample(1.,1.,2.);
    auto results=sim.collect_results(&reducer);
    bool ok=results.empty();
    if (!rank) {
      auto const& density=results.at("density_0");
      auto const& large=results.at("g2w_re_0_0");
      auto count=active<0 ? 0 : 3;
      ok=density.size()==1 && large.size()==4 && density.count()==count && large.count()==count;
      if (active>=0) ok=ok && density.mean<double>()(0)==2. && density.stderror<double>()(0)==0.
        && large.mean<double>()(0)==2. && large.stderror<double>()(0)==0.;
    }
    require(ok, "empty CT-HYB replica changed a mixed-family root result");
  }
  fixture unequal(run,rank);
  for (int i=0; i<(rank ? 5 : 3); ++i)
    unequal.sample(rank && i>=3 ? -1. : 1., rank && i>=4 ? -1. : 1., rank ? 5. : 2.);
  auto results=unequal.collect_results(&reducer);
  bool ok=results.empty();
  if (!rank) {
    auto const& density=results.at("density_0");
    auto const& large=results.at("g2w_re_0_0");
    // Small joint batches: pseudovalues 1 (three), 8 (three), -.4 (two).
    // Large joint moments: ratio 21/6; residuals +/-1.5 give error^2=4/7.
    ok=results.at("Sign").count()==8 && results.at("Sign").mean<double>()(0)==.5
      && density.count()==8 && std::abs(density.mean<double>()(0)-3.275)<1.e-12
      && std::abs(density.stderror<double>()(0)-std::sqrt(1.955625))<1.e-12
      && large.count()==8 && std::abs(large.mean<double>()(0)-3.5)<1.e-12
      && std::abs(large.mean<double>()(3)-14.)<1.e-12
      && std::abs(large.stderror<double>()(0)-std::sqrt(4./7.))<1.e-12;
  }
  require(ok, "CT-HYB reduced corrected replicas or used the wrong measurement-cadence denominator");
  fixture zero_sign(run,rank);
  zero_sign.sample(rank ? -1. : 1., rank ? -1. : 1.,2.);
  rejects([&] { zero_sign.collect_results(&reducer); });
  require(zero_sign.unchanged(rank ? -1. : 1.), "failed CT-HYB analysis changed live measurements");
  for (int mismatch=0; mismatch<4; ++mismatch) {
    fixture invalid(run,rank);
    if (rank) {
      if (!mismatch) invalid.change_name();
      else if (mismatch==1) invalid.change_shape();
      else if (mismatch==2) invalid.change_family();
      else invalid.change_sign();
    }
    rejects([&] { invalid.collect_results(&reducer); });
  }
}
void driver(alps::run_configuration run, int rank) {
  auto const path=run.output["results"].as<std::string>();
  auto const text=std::filesystem::path("cthyb-mpi-text");
  require(rank || (!std::filesystem::exists(path) && !std::filesystem::exists(text)),
          "CT-HYB test output already exists");
  alps::solvers::cthyb(run);
  bool ok=true;
  std::string previous;
  if (!rank) try {
    alps::hdf5::archive archive(path);
    alps::alea::hdf5_serializer serializer(archive,"/impurity/simulation/results");
    alps::alea::batch_result<double> sign;
    alps::alea::var_result<double> large;
    deserialize(serializer,"Sign",sign);
    deserialize(serializer,"g2w_re_0_0",large);
    std::vector<double> green;
    archive["/impurity/G_tau/0/mean/value"]>>green;
    std::string application;
    archive["/run_config/application"]>>application;
    ok=sign.count()>0 && sign.mean()(0)==1. && large.count()==sign.count() && large.size()==4
      && application=="cthyb" && green.size()==3 && green.front()+green.back()==-1.;
    archive.close();
    std::ifstream file(path,std::ios::binary);
    previous.assign(std::istreambuf_iterator<char>(file),{});
    std::filesystem::create_directories(text/"orders.dat");
  } catch (std::exception const& error) { std::cerr<<error.what()<<'\n'; ok=false; }
  require(ok, "native CT-HYB driver did not publish both result families, Green functions and provenance");
  // Valid directory, but text publication fails after native result serialization.
  run.output["text"]=true;
  run.output["text_directory"]=std::filesystem::absolute(text).string();
  rejects([&] { alps::solvers::cthyb(run); });
  if (!rank) {
    std::ifstream file(path,std::ios::binary);
    ok=previous==std::string(std::istreambuf_iterator<char>(file),{});
    std::filesystem::remove(path);
    std::filesystem::remove_all(text);
  }
  require(ok, "failed CT-HYB publication changed the previous HDF5 output");
}
}
int main(int argc, char** argv) {
  MPI_Init(&argc,&argv);
  int rank,size,status=0;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  try {
    require(size==2, "run CT-HYB contract with two ranks");
    require(rank || !std::filesystem::exists("cthyb-mpi-delta.dat"), "CT-HYB test input already exists");
    if (!rank) std::ofstream("cthyb-mpi-delta.dat")<<"0 -0.5 -0.5\n1 -0.5 -0.5\n2 -0.5 -0.5\n";
    MPI_Barrier(MPI_COMM_WORLD);
    auto run=configuration();
    collection(run,MPI_COMM_WORLD);
    MPI_Comm reversed;
    MPI_Comm_split(MPI_COMM_WORLD,0,1-rank,&reversed);
    collection(run,reversed);
    MPI_Comm_free(&reversed);
    driver(run,rank);
    if (!rank) std::filesystem::remove("cthyb-mpi-delta.dat");
  } catch (std::exception const& error) {
    if (!rank) std::cerr<<error.what()<<'\n';
    status=1;
  }
  MPI_Finalize();
  return status;
}
