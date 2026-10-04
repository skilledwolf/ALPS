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

#include "run_config.h"
#include "dmft_schema.hpp"
#include <alps/run_config.hpp>
#include <alps/utility/copyright.hpp>
#include <iostream>
#include <string>
#include <utility>
#ifdef ALPS_HAVE_MPI
#include <boost/mpi/environment.hpp>
#endif

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
#ifdef ALPS_HAVE_MPI
    boost::mpi::environment environment(argc, argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
    const int rank=0;
#endif
    alps::run_configuration run;
    std::string failure;
    try {
      run=alps::load_run_configuration(filename,alps::dmft::hirschfye_schema);
      if (validate) alps::dmft::prepare_hirschfye_run(run);
    } catch (std::exception const& error) { failure=error.what(); }
    alps::dmft::agree_hirschfye_failure(failure);
    if (validate) {
      if (!rank) std::cout << "Valid Hirsch-Fye configuration: " << filename << '\n';
      return 0;
    }
    if (!rank) {
      std::cout << "ALPS Hirsch-Fye solver for the single site impurity problem.\n\n";
      alps::print_copyright(std::cout);
      std::cout << "****************************************************************\n"
                 "* Recommended citation in scientific publications:             *\n"
                 "* We used the ALPS [1] implementation [2] of the Hirsch-Fye    *\n"
                 "* [3] impurity solver.                                         *\n"
                 "* [1] JSTAT (2011) P05001; [2] CPC 182, 1078 (2011); [3] PRL   *\n"
                 "* 56, 2521 (1986).                                             *\n"
                 "****************************************************************\n";
    }
    alps::dmft::run_hirschfye(std::move(run));
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "hirschfye: " << error.what() << '\n';
    return 1;
  }
}
