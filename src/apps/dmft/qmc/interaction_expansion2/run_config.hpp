// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#pragma once
#include <alps/run_config.hpp>
#include "../green_function.h"
// Shared input decoding for the simulation and result postprocessing.
void read_ctint_bare_green(const alps::params &parameters, const alps::params &input,
                          matsubara_green_function_t &green);
namespace alps::ctint {
std::string schema_for_run(const std::filesystem::path &file);
void validate_interaction(const params &parameters, const params &input);
void validate_execution(const params &execution);
void agree_failure(const std::string &error);
}
