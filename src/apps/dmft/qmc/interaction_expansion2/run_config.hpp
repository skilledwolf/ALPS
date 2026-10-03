// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#pragma once
#include <alps/run_config.hpp>
#include "green_function.h"
// Shared input decoding for the simulation and result postprocessing.
void read_ctint_bare_green(const alps::params &parameters, const alps::params &input,
                          matsubara_green_function_t &green);
