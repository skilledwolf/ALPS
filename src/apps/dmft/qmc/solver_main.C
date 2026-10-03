/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2010 by Emanuel Gull <gull@phys.columbia.edu>
 *                              Philipp Werner <werner@itp.phys.ethz.ch>,
 *                              Sebastian Fuchs <fuchs@theorie.physik.uni-goettingen.de>
 *                              Matthias Troyer <troyer@comp-phys.org>
 *
 *
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
 *
 *****************************************************************************/

/* $Id: solver_main.C 285 2008-01-03 15:34:39Z gullc $ */

#include "hirschfyesim.h"
#include "dmft_schema.hpp"
#include <alps/hdf5/complex.hpp>
#include <alps/run_config.hpp>
#include <alps/utility/copyright.hpp>
#include <iostream>
#include <string>

int main(int argc, char** argv) {
  try {
    bool validate = false, show_schema = false;
    std::string filename;
    for (int i = 1; i < argc; ++i) {
      const std::string argument(argv[i]);
      if (argument == "--help" || argument == "-h") {
        std::cout << "Usage: hirschfye [--validate] run.toml | hirschfye --schema\n";
        return 0;
      }
      if (argument == "--schema") show_schema = true;
      else if (argument == "--validate") validate = true;
      else if (argument.empty() || argument.front() == '-')
        throw std::invalid_argument("Unknown option: " + argument);
      else if (filename.empty()) filename = argument;
      else throw std::invalid_argument("Expected one TOML run file");
    }
    if (show_schema) { std::cout << alps::dmft::hirschfye_schema; return 0; }
    if (filename.empty()) throw std::invalid_argument("No TOML run file specified");
    auto run = alps::load_run_configuration(filename, alps::dmft::hirschfye_schema);
    if (run.parameters["BETA"].as<double>() <= 0.)
      throw std::invalid_argument("Hirsch-Fye BETA must be positive");
    matsubara_green_function_t g0(run.parameters["NMATSUBARA"].as<unsigned int>(),
                                  run.parameters["SITES"].as<unsigned int>(),
                                  run.parameters["FLAVORS"].as<unsigned int>());
    {
      alps::hdf5::archive input(run.input["g0"].as<std::string>(), "r");
      read_flavor_vectors(input, "/G0", g0);
    }
    if (validate) {
      std::cout << "Valid Hirsch-Fye configuration: " << filename << '\n';
      return 0;
    }
    std::cout << "ALPS Hirsch-Fye solver for the single site impurity problem.\n\n";
    alps::print_copyright(std::cout);
    std::cout << "****************************************************************\n"
                 "* Recommended citation in scientific publications:             *\n"
                 "* We used the ALPS [1] implementation [2] of the Hirsch-Fye    *\n"
                 "* [3] impurity solver.                                         *\n"
                 "* [1] JSTAT (2011) P05001; [2] CPC 182, 1078 (2011); [3] PRL   *\n"
                 "* 56, 2521 (1986).                                             *\n"
                 "****************************************************************\n";
    alps::scheduler::BasicFactory<HirschFyeSim,HirschFyeRun> factory;
    alps::ImpuritySolver solver(factory, run, argc, argv);
    solver.solve_omega(g0, run.parameters);
    alps::hdf5::archive output(run.output["results"].as<std::string>(), "a");
    output["/run_config"] << run;
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "hirschfye: " << error.what() << '\n';
    return 1;
  }
}
