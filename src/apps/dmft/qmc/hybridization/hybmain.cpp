// Copyright (C) 2012 Emanuel Gull, Hartmut Hafermann.
// Modifications (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#include <alps/cthyb.hpp>
#include <alps/solvers.hpp>
#include <iostream>
#ifdef ALPS_HAVE_MPI
#include <boost/mpi/environment.hpp>
#endif

int main(int argc, char** argv) {
  try {
    bool validate=false;
    std::string filename;
    for(int i=1;i<argc;++i){
      const std::string arg=argv[i];
      if(arg=="--schema") { std::cout << alps::cthyb::schema(); return 0; }
      if(arg=="--help" || arg=="-h"){
        std::cout<<"Usage: hybridization [--validate] run.toml | hybridization --schema\n"
                 <<"Input, output, and execution settings belong in the TOML run file.\n";
        return 0;
      }
      if(arg=="--validate") validate=true;
      else if(arg.empty() || arg.front()=='-') throw std::invalid_argument("Unknown option: "+arg);
      else if(filename.empty()) filename=arg;
      else throw std::invalid_argument("Expected one TOML run file");
    }
    if(filename.empty()) throw std::invalid_argument("No TOML run file specified");
    auto run=alps::load_run_configuration(filename,alps::cthyb::schema());
    alps::cthyb::prepare_run(run);
    if(validate){std::cout<<"Valid CT-HYB configuration: "<<filename<<'\n'; return 0;}
#ifdef ALPS_HAVE_MPI
    boost::mpi::environment env(argc, argv);
#endif
    alps::solvers::cthyb(run);
    return 0;
  } catch(const std::exception& error) {
    std::cerr<<"hybridization: "<<error.what()<<'\n';
    return 1;
  }
}
