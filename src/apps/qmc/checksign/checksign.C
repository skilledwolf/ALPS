/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2003-2006 by Matthias Troyer <troyer@itp.phys.ethz.ch>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include "lattice_model.hpp"
#include <alps/model.h>
#include <alps/lattice.h>
#include <alps/utility/copyright.hpp>
#include <cstring>
#include <iostream>
#include <string>

// Each argument is a TOML run file with the lattice and model in [parameters]
// and their libraries in input.lattice_library and input.model_library.
constexpr char run_schema[] = R"toml(application = "checksign"
schema_version = 1
[parameters.LATTICE]
type = "string"
[parameters.GRAPH]
type = "string"
[parameters.MODEL]
type = "string"
required = true
[input.lattice_library]
type = "path"
[input.model_library]
type = "path"
[output]
[execution]
)toml";

int main(int argc, char** argv)
{
  try {
    std::cout << "ALPS application to check for a sign problem in a quantum model\n"
              << "  available from http://alps.comp-phys.org/\n"
              << "  copyright (c) 2003-2007 by Matthias Troyer <troyer@comp-phys.org>\n\n";
    alps::print_copyright(std::cout);

    if (argc<2) {
      std::cerr << "Usage: " << argv[0] << " [-l] run.toml [run.toml ...]\n";
      return 1;
    }
    for (int i=1;i<argc;++i) {
      if (!std::strcmp(argv[i],"-l")) {
        alps::print_license(std::cout);
        continue;
      }
      auto run = alps::load_run_configuration(argv[i], lattice_model::schema(argv[i], run_schema));
      lattice_model::resolve_libraries(run);
      const alps::Parameters p = lattice_model::parameters(run);
      alps::graph_helper<> lattice(p);
      alps::model_helper<> models(lattice, p);
      std::cout << argv[i] << (alps::has_sign_problem(models.model(),lattice,p) ? ": SIGN PROBLEM\n" : ": OK\n");
    }
  }
  catch (std::exception& exc) {
    std::cerr << exc.what() << "\n";
    return 1;
  }
}
