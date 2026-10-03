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
#include <alps/run_config.hpp>
#include <iostream>
#ifdef ALPS_HAVE_MPI
#include <boost/mpi/environment.hpp>
#endif

int main(int argc, char** argv) {
  try {
    bool validate = false;
    std::string file;
    for (int i = 1; i < argc; ++i) {
      const std::string argument = argv[i];
      if (argument == "--schema") { std::cout << alps::ctint::schema(); return 0; }
      if (argument == "--help" || argument == "-h") {
        std::cout << "Usage: interaction [--validate] run.toml | interaction --schema\n"
                  << "Use [parameters], [input], [output], and [execution] TOML sections.\n";
        return 0;
      }
    }
    for (int i = 1; i < argc; ++i) {
      const std::string argument = argv[i];
      if (argument == "--validate") validate = true;
      else if (argument.empty() || argument[0] == '-')
        throw std::invalid_argument("unknown option: " + argument);
      else if (file.empty()) file = argument;
      else throw std::invalid_argument("expected one TOML run file");
    }
    if (file.empty()) throw std::invalid_argument("No TOML run file specified");
#ifdef ALPS_HAVE_MPI
    boost::mpi::environment environment(argc, argv);
#endif
    auto run = alps::load_run_configuration(file, alps::ctint::schema());
    if (validate) {
      alps::ctint::prepare_run(run);
      std::cout << "Valid CT-INT configuration: " << file << '\n';
    } else alps::solvers::ctint(run);
    return 0;
  } catch (const std::exception &error) {
    std::cerr << "interaction: " << error.what() << '\n';
    return 1;
  }
}
