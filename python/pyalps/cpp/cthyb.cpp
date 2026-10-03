// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
#include "run_config.hpp"
#include <alps/cthyb.hpp>

NB_MODULE(cthyb, module) {
  pyalps::bind_configured_solver(module, alps::cthyb::schema(), alps::cthyb::prepare_run,
                                 alps::solvers::cthyb);
}
