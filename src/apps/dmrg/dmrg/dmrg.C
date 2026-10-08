/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2006 -2010 by Adrian Feiguin <afeiguin@uwyo.edu>
*                             Matthias Troyer <troyer@itp.phys.ethz.ch>
* Modifications (C) 2026 ALPS Collaboration
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include "dmrg.h"
#include "lattice_model.hpp"
#include "schema.hpp"
#include <complex>
#include <cstdint>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
struct application {
  static void print_copyright(std::ostream& out) { print_dmrg_copyright(out); }

  template <class F>
  static void with_task(alps::run_configuration const& run, F const& f)
  {
    alps::SymbolTable p = lattice_model::parameters(run);
    // The sweep schedule is read as a comma-separated list.
    if (run.parameters.exists("STATES")) {
      std::string states;
      for (std::int64_t n : run.parameters["STATES"].as<std::vector<std::int64_t> >())
        states += (states.empty() ? "" : ",") + std::to_string(n);
      p["STATES"] = states;
    }
    const std::string directory = run.execution.exists("temporary_directory")
      ? run.execution["temporary_directory"].as<std::string>() : alps::temp_directory_path().string();
    if (!std::filesystem::is_directory(directory))
      throw std::invalid_argument("execution.temporary_directory must name an existing directory");
    if (p.value_or_default("COMPLEX", false)) {
      DMRGTask<std::complex<double> > task(p, directory);
      f(task);
    } else {
      DMRGTask<double> task(p, directory);
      f(task);
    }
  }
};
} // namespace

int main(int argc, char** argv)
{
  return lattice_model::main<application>(argc, argv, "dmrg", run_schema);
}
