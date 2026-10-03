// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#include <alps/cthyb.hpp>
#include <alps/solvers.hpp>
#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/vector.hpp>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#ifdef ALPS_HAVE_MPI
#include <boost/mpi/environment.hpp>
#endif

namespace {
void check(bool ok,const char* message) { if(!ok) throw std::runtime_error(message); }
template<class F> void rejects(F f,const char* expected) {
  try { f(); } catch(const std::exception& e) {
    if(std::string(e.what()).find(expected)!=std::string::npos) return;
    throw std::runtime_error(std::string("Wrong rejection: ")+e.what());
  }
  throw std::runtime_error(std::string("Accepted invalid configuration: ")+expected);
}
alps::params minimal() {
  alps::params p;
  p["BETA"]=2.; p["N_ORBITALS"]=2; p["N_TAU"]=2; p["N_MEAS"]=1;
  p["THERMALIZATION"]=0; p["SWEEPS"]=10; p["U"]=1.;
  return p;
}
}
int main(int argc,char** argv) {
#ifdef ALPS_HAVE_MPI
  boost::mpi::environment environment(argc,argv);
#endif
  const auto directory=std::filesystem::temp_directory_path()/
    ("alps-cthyb-config-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
  std::filesystem::create_directory(directory);
  struct cleanup { std::filesystem::path path; ~cleanup(){std::error_code error;std::filesystem::remove_all(path,error);} } guard{directory};
  try {
    const auto p=alps::cthyb::prepare_parameters(minimal());
    check(p["MEASURE_time"].as<bool>() && !p["MEASURE_freq"].as<bool>(),"Measurement defaults changed");
    auto bad=minimal(); bad["N_MEAS"]=0;
    rejects([&]{alps::cthyb::prepare_parameters(bad);},"N_MEAS");
    bad=minimal(); bad["BETA"]=0.;
    rejects([&]{alps::cthyb::prepare_parameters(bad);},"BETA");
    bad=minimal(); bad["MEASURE_freq"]=true;
    rejects([&]{alps::cthyb::prepare_parameters(bad);},"N_MATSUBARA");
    bad=minimal(); bad["MEASURE_g2w"]=true; bad["N_w2"]=3; bad["N_W"]=1;
    rejects([&]{alps::cthyb::prepare_parameters(bad);},"even");
    bad["N_w2"]=4; bad["COMPUTE_VERTEX"]=true;
    rejects([&]{alps::cthyb::prepare_parameters(bad);},"MEASURE_freq");
    bad["MEASURE_freq"]=true; bad["N_MATSUBARA"]=1;
    rejects([&]{alps::cthyb::prepare_parameters(bad);},"N_MATSUBARA >=");
    bad["N_MATSUBARA"]=2;
    alps::cthyb::prepare_parameters(bad);
    bad=minimal(); bad["MAX_TIME"]=10;
    rejects([&]{alps::cthyb::prepare_parameters(bad);},"unknown key");
    bad=minimal(); bad["N_ORBITALS"]=50000;
    rejects([&]{alps::cthyb::prepare_parameters(bad);},"index range");
    const std::string delta_text="0 -0.5 -0.5\n1 -0.5 -0.5\n2 -0.5 -0.5\n";
    std::ofstream(directory/"delta.dat")<<delta_text;
    {
      struct restore_directory {
        std::filesystem::path previous=std::filesystem::current_path();
        ~restore_directory(){std::error_code error;std::filesystem::current_path(previous,error);}
      } restore;
      std::filesystem::current_path(directory);
      const auto base=std::filesystem::current_path();
      const auto input_unchanged=[&]{
        std::ifstream input(base/"delta.dat");
        check(std::string(std::istreambuf_iterator<char>(input),{})==delta_text,
              "Rejected path collision changed the scientific input");
      };
      alps::run_configuration relative;
      relative.parameters=minimal(); relative.input["delta"]="delta.dat";
      relative.output["results"]="delta.dat";
      rejects([&]{alps::cthyb::prepare_run(relative);},"must not replace input.delta");
      input_unchanged();
      std::error_code symlink_error;
      std::filesystem::create_symlink("delta.dat","delta-alias.dat",symlink_error);
      if(!symlink_error){
        relative.output["results"]="delta-alias.dat";
        rejects([&]{alps::cthyb::prepare_run(relative);},"must not replace input.delta");
        input_unchanged();
      } else {
#ifdef _WIN32
        std::cout<<"Symlink collision check unavailable: "<<symlink_error.message()<<'\n';
#else
        throw std::filesystem::filesystem_error("Create input alias",base/"delta-alias.dat",symlink_error);
#endif
      }
      relative.output["results"]="relative-result.h5";
      alps::cthyb::prepare_run(relative);
      check(relative.input["delta"].as<std::string>()==(base/"delta.dat").string(),
            "Prepared scientific input must be anchored to the working directory");
      check(relative.output["results"].as<std::string>()==(base/"relative-result.h5").string(),
            "Prepared results must be anchored to the working directory");
      std::filesystem::create_directory(base/"other");
      std::filesystem::current_path(base/"other");
      alps::cthyb::prepare_run(relative);
      check(relative.input["delta"].as<std::string>()==(base/"delta.dat").string() &&
            relative.output["results"].as<std::string>()==(base/"relative-result.h5").string(),
            "A working-directory change reinterpreted prepared paths");
      input_unchanged();
    }
    alps::run_configuration run;
    run.parameters=minimal(); run.input["delta"]=(directory/"delta.dat").string();
    run.output["results"]=(directory/"result.h5").string();
    alps::cthyb::prepare_run(run);
    check(run.application=="cthyb" && run.schema_version==1,"Application identity must come from schema");
    check(run.execution["time_limit"].as<int>()==0,"Finite-sweep run must not require a time limit");
    auto invalid=run; invalid.output["results"]=invalid.input["delta"];
    rejects([&]{alps::cthyb::prepare_run(invalid);},"must not replace input.");
    invalid=run; invalid.output["text"]=true;
    invalid.output["text_directory"]=directory.string();
    invalid.output["results"]=(directory/"simulation.dat").string();
    rejects([&]{alps::cthyb::prepare_run(invalid);},"replace HDF5 results");
    std::ofstream(directory/"positive.dat")<<"0 0.5 -0.5\n1 -0.5 -0.5\n2 -0.5 -0.5\n";
    invalid=run; invalid.input["delta"]=(directory/"positive.dat").string();
    rejects([&]{alps::cthyb::prepare_run(invalid);},"nonpositive");
    std::ofstream(directory/"short.dat")<<"0 -0.5 -0.5\n";
    invalid.input["delta"]=(directory/"short.dat").string();
    rejects([&]{alps::cthyb::prepare_run(invalid);},"malformed");
    std::ofstream(directory/"coordinate.dat")<<"0 -0.5 -0.5\n0.5 -0.5 -0.5\n2 -0.5 -0.5\n";
    invalid.input["delta"]=(directory/"coordinate.dat").string();
    rejects([&]{alps::cthyb::prepare_run(invalid);},"coordinate grid");
    {
      alps::hdf5::archive ar((directory/"delta.h5").string(),"w");
      ar["/Delta_0"]<<std::vector<double>{-0.5,-0.5,-0.5};
      ar["/Delta_1"]<<std::vector<double>{-0.5,-0.5};
    }
    invalid=run; invalid.input["delta"]=(directory/"delta.h5").string(); invalid.input["delta_format"]="hdf5";
    rejects([&]{alps::cthyb::prepare_run(invalid);},"length 3");
    {
      alps::hdf5::archive ar((directory/"delta.h5").string(),"a");
      ar["/Delta_1"]<<std::vector<double>{-0.5,std::numeric_limits<double>::infinity(),-0.5};
    }
    rejects([&]{alps::cthyb::prepare_run(invalid);},"finite");
    std::ofstream(directory/"mu.dat")<<"0.2 0.3";
    run.input["chemical_potential"]=(directory/"mu.dat").string();
    alps::cthyb::prepare_run(run);
    run.parameters["MEASURE_time"]=false;
    alps::solvers::cthyb(run);
    {
      alps::hdf5::archive ar(run.output["results"].as<std::string>(),"r");
      std::string application;
      ar["/run_config/application"]>>application;
      check(application=="cthyb","Results lost application identity");
      alps::params execution,parameters;
      ar["/run_config/execution"]>>execution;
      ar["/parameters"]>>parameters;
      check(execution["time_limit"].as<int>()==0,"Finite-sweep execution metadata changed");
      check(!parameters.exists("DELTA") && !parameters.exists("MAX_TIME"),"Orchestration leaked into scientific parameters");
      check(ar.is_group("/simulation/results/Sign"),"Simulation did not produce measurements");
    }
    std::cout<<"CT-HYB run configuration contracts passed\n";
    return 0;
  } catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
