// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
#ifndef ALPS_SOLVERS_HPP
#define ALPS_SOLVERS_HPP

#include <string>

namespace alps {
namespace params_ns { class dictionary; }
using params = params_ns::dictionary;
struct run_configuration;

// Link the corresponding ALPS::maxent, ALPS::cthyb, or ALPS::ctint target.
// Solvers run synchronously, write to output_file, and propagate exceptions.
// With an MPI-enabled SDK, the caller must initialize MPI before CT-QMC calls;
// all ranks in MPI_COMM_WORLD must participate.
namespace maxent { struct data; }
namespace solvers {
bool maxent(params const& parameters, const alps::maxent::data& data,
            std::string const& output_file, int time_limit=60, bool text_output=false);
bool maxent(run_configuration const& run);
void cthyb(run_configuration const& run);
void ctint(run_configuration const& run);
}
}

#endif
