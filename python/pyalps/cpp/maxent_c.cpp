// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#include "run_config.hpp"
#include <alps/maxent.hpp>
NB_MODULE(maxent_c,module) {
    pyalps::bind_configured_solver(module, alps::maxent::schema(), alps::maxent::prepare_run,
        static_cast<bool (*)(const alps::run_configuration &)>(&alps::solvers::maxent));
}
