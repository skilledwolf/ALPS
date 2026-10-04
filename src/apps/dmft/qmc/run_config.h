// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#pragma once
#include <alps/run_config.hpp>
#include "green_function.h"
#include <string>
namespace alps::dmft {
// External executables (in the ALPS bin directory, or a custom path) or the
// in-process scheduler CT-INT. A custom executable supplies input.solver_schema
// and execution.solver_input.
enum class solver_kind { hybridization, interaction, hirsch_fye, interaction_expansion, custom };
solver_kind selected_solver(const run_configuration& run);
bool external(solver_kind kind);
// True when the solver receives a hybridization function rather than G0.
bool receives_delta(const run_configuration& run);
// Schema text of the selected solver; a relative custom schema resolves against base.
std::string solver_schema(const run_configuration& run, const std::filesystem::path& base = {});
// DMFT parameters under the solver's names, restricted to the keys it declares.
params solver_parameters(const params& dmft, std::string_view schema);
// Scientific inputs, such as a retarded interaction, passed unchanged to a solver
// that declares them.
params solver_inputs(const params& dmft, std::string_view schema);
std::string schema_for_run(const std::filesystem::path& file = {});
run_configuration load_run(const std::filesystem::path& file);
params prepare_hirschfye_parameters(const params&);
// Resolves and validates the run, returning its once-read scientific input.
matsubara_green_function_t prepare_hirschfye_run(run_configuration&);
// Internal CLI/driver failure boundary; collective in an MPI build.
void agree_hirschfye_failure(const std::string&);
void run_hirschfye(run_configuration);
}
