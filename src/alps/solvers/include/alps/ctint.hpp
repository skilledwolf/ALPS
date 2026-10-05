// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#pragma once
#include <alps/run_config.hpp>
#include <alps/solvers.hpp>
namespace alps::ctint {
std::string schema(const params &parameters = {});
params prepare_parameters(const params &supplied);
void prepare_run(run_configuration &run);
}
