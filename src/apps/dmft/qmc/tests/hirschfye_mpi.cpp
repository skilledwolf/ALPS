// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#include "hirschfyesim.h"
#include "run_config.h"
#include <alps/alea/mpi.hpp>
#include <alps/alea/hdf5.hpp>
#include <alps/hdf5/complex.hpp>
#include <alps/hdf5/vector.hpp>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>

namespace {
void require(bool value, char const* message) {
  int agreed=value;
  MPI_Allreduce(MPI_IN_PLACE,&agreed,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  if (!agreed) throw std::runtime_error(message);
}
template<class F> void rejects(F operation) {
  bool failed=false;
  try { operation(); } catch (std::exception const&) { failed=true; }
  require(failed,"invalid Hirsch-Fye operation was not rejected on every rank");
}
alps::run_configuration configuration() {
  alps::run_configuration run;
  run.parameters["BETA"]=2.; run.parameters["U"]=0.;
  run.parameters["N"]=4; run.parameters["NMATSUBARA"]=4;
  run.parameters["EPSSQ_0"]=0.; run.parameters["EPSSQ_1"]=0.;
  // Finish at the first collective while retaining the default check schedule.
  run.parameters["SWEEPS"]=1; run.parameters["THERMALIZATION"]=0;
  run.input["g0"]="hirschfye-mpi-input.h5";
  run.output["results"]="hirschfye-mpi-result.h5";
  run.execution["bins"]=16;
  return run;
}
struct fixture : HirschFyeRun {
  using HirschFyeRun::HirschFyeRun;
  void sample(double sign,double value) {
    measurements[0]<<sign;
    alps::alea::column<double> green(6);
    green.head(4).setConstant(-sign*value);
    green(4)=sign*(value-1.); green(5)=sign;
    measurements[1]<<green; measurements[2]<<green;
  }
  void change_shape() { measurements[1].set_size(7); }
  bool unchanged(double sign) const {
    auto raw=measurements[1].result();
    return raw.count()==1 && raw.mean()(0)==-2*sign && raw.mean()(5)==sign;
  }
};
void collection(alps::run_configuration const& run,matsubara_green_function_t const& g0,MPI_Comm comm) {
  int rank;
  MPI_Comm_rank(comm,&rank);
  alps::alea::mpi_reducer reducer(comm);
  for (int active : {0,1,-1}) {
    fixture simulation(run,g0,rank);
    if (rank==active) for (int i=0;i<3;++i) simulation.sample(1.,.5);
    auto results=simulation.collect_results(&reducer);
    bool ok=results.empty();
    if (!rank) {
      auto const& green=results.at("G_meas_up");
      ok=green.count()==(active<0 ? 0 : 3) && green.size()==5;
      if (active>=0) ok=ok && green.mean()(0)==-.5 && green.mean()(4)==-.5
        && green.stderror()(0)==0.;
    }
    require(ok,"empty Hirsch-Fye replica changed the pooled physical result");
  }
  fixture unequal(run,g0,rank);
  for (int i=0;i<(rank ? 5 : 3);++i)
    unequal.sample(rank && i>=3 ? -1. : 1.,rank ? 5. : 2.);
  auto results=unequal.collect_results(&reducer);
  bool ok=results.empty();
  if (!rank) {
    auto const& green=results.at("G_meas_up");
    // Independent partial-bin pseudovalues: -1 (three), -8 (three), .4 (two).
    ok=results.at("Sign").count()==8 && results.at("Sign").mean()(0)==.5
      && green.count()==8 && std::abs(green.mean()(0)+3.275)<1.e-12
      && std::abs(green.stderror()(0)-std::sqrt(1.955625))<1.e-12
      && std::abs(green.mean()(0)+green.mean()(4)+1.)<1.e-12;
  }
  require(ok,"Hirsch-Fye corrected replicas before pooling raw signed batches");
  fixture zero_sign(run,g0,rank);
  zero_sign.sample(rank ? -1. : 1.,2.);
  rejects([&] { zero_sign.collect_results(&reducer); });
  require(zero_sign.unchanged(rank ? -1. : 1.),"failed Hirsch-Fye analysis changed live measurements");
  fixture invalid(run,g0,rank);
  if (rank) invalid.change_shape();
  rejects([&] { invalid.collect_results(&reducer); });
}
void driver(alps::run_configuration run,int rank) {
  const auto directory=std::filesystem::absolute("hirschfye-mpi-output");
  const auto path=directory/"result.h5";
  bool ok=true;
  if (!rank) {
    ok=!std::filesystem::exists(directory);
    if (ok) std::filesystem::create_directory(directory);
  }
  require(ok,"Hirsch-Fye test output already exists");
  run.output["results"]=path.string();
  alps::dmft::run_hirschfye(run);
  std::string previous;
  HirschFyeRun::results_type snapshots;
  if (!rank) try {
    alps::hdf5::archive archive(path.string());
    alps::alea::hdf5_serializer serializer(archive,"/simulation/results");
    alps::alea::batch_result<double> sign,green;
    deserialize(serializer,"Sign",sign); deserialize(serializer,"G_meas_up",green);
    std::vector<double> tau;
    std::vector<std::complex<double>> omega;
    archive["/G_tau/0/mean/value"]>>tau;
    archive["/G_omega/0/mean/value"]>>omega;
    std::string application;
    archive["/run_config/application"]>>application;
    ok=application=="hirschfye" && sign.count()==2 && sign.mean()(0)==1.
      && green.count()==2 && green.size()==5 && tau.size()==5 && omega.size()==4;
    for (auto value : tau) ok=ok && std::abs(value+.5)<1.e-12;
    for (std::size_t i=0;i<omega.size();++i)
      ok=ok && std::abs(omega[i]-std::complex<double>(0.,-2./((2*i+1)*std::acos(-1.))))<1.e-12;
    for (auto const* name : {"Sign","G_meas_up","G_meas_down"}) {
      alps::alea::batch_result<double> result;
      deserialize(serializer,name,result);
      snapshots.emplace(name,std::move(result));
    }
    archive.close();
    std::ifstream file(path,std::ios::binary);
    previous.assign(std::istreambuf_iterator<char>(file),{});
  } catch (std::exception const& error) { std::cerr<<error.what()<<'\n'; ok=false; }
  require(ok,"native Hirsch-Fye driver did not publish canonical batches, exact free Green functions and provenance");
  auto unavailable=run;
  if (rank) unavailable.input["g0"]="hirschfye-mpi-missing-input.h5";
  rejects([&] { alps::dmft::run_hirschfye(unavailable); });
  if (!rank) {
    std::ifstream file(path,std::ios::binary);
    ok=previous==std::string(std::istreambuf_iterator<char>(file),{});
  }
  require(ok,"rank-local Hirsch-Fye initialization failure changed the prior archive");
  bool serialized=false,failed=false;
  if (!rank) {
    // Use real production snapshots and the same checked publication boundary.
    // A missing spin makes scientific evaluation fail after all batches write.
    try {
      alps::hdf5::save_checkpoint(path.string(),[&](alps::hdf5::archive& archive) {
        archive["/parameters"]<<run.parameters;
        alps::alea::hdf5_serializer serializer(archive,"/simulation/results");
        for (auto const& [name,result] : snapshots) serialize(serializer,name,result);
        serialized=true;
        auto incomplete=snapshots;
        incomplete.erase("G_meas_down");
        HirschFyeRun::get_result(incomplete,run.parameters);
      });
    } catch (std::exception const&) { failed=true; }
    std::ifstream file(path,std::ios::binary);
    ok=previous==std::string(std::istreambuf_iterator<char>(file),{});
  }
  require(rank || (serialized && failed),"Hirsch-Fye scientific evaluation did not fail after serialization");
  require(ok,"failed Hirsch-Fye publication changed the prior archive");
  if (!rank) std::filesystem::remove_all(directory);
}
}
int main(int argc,char** argv) {
  MPI_Init(&argc,&argv);
  int rank,size,status=0;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank); MPI_Comm_size(MPI_COMM_WORLD,&size);
  try {
    require(size==2,"run Hirsch-Fye contract with two ranks");
    const auto input=std::filesystem::path("hirschfye-mpi-input.h5");
    require(rank || !std::filesystem::exists(input),"Hirsch-Fye test input already exists");
    if (!rank) {
      alps::hdf5::archive archive(input.string(),"w");
      std::vector<std::complex<double>> values(4);
      for (std::size_t i=0;i<values.size();++i) values[i]={0.,-2./((2*i+1)*std::acos(-1.))};
      archive["/G0_0"]<<values; archive["/G0_1"]<<values;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    auto run=configuration();
    const auto g0=alps::dmft::prepare_hirschfye_run(run);
    collection(run,g0,MPI_COMM_WORLD);
    MPI_Comm reversed;
    MPI_Comm_split(MPI_COMM_WORLD,0,1-rank,&reversed);
    collection(run,g0,reversed);
    MPI_Comm_free(&reversed);
    driver(run,rank);
    if (!rank) std::filesystem::remove(input);
  } catch (std::exception const& error) {
    if (!rank) std::cerr<<error.what()<<'\n';
    status=1;
  }
  MPI_Finalize();
  return status;
}
