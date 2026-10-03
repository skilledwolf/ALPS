// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
#include "run_config.hpp"
#include <alps/ctint.hpp>

NB_MODULE(ctint, module) {
  pyalps::bind_configured_solver(module, alps::ctint::schema(), alps::ctint::prepare_run,
                                 alps::solvers::ctint);
}
