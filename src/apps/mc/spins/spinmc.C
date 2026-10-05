// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#include "spinsim.h"
#include "schema.hpp"
#include "../native_driver.hpp"

int main(int argc, char** argv) {
    return native_mc::main<spinmc::simulation>(argc, argv, "spinmc", spinmc_schema,
        [](std::string const& key, toml::node const&) -> char const* {
            if (key == "J'" || key == "D'")
                throw std::invalid_argument("Use J/D or numbered J<bond>/D<site> matrix arrays");
            if (key.size() > 1 && key.find_first_not_of("0123456789", 1) == std::string::npos) {
                if (key[0] == 'J' || key[0] == 'D') return "float64[]";
                if (key[0] == 'S') return "float64";
            }
            return nullptr;
        },
        [](alps::params& p, alps::run_configuration const& run) {
            if (p.exists("T") == p.exists("beta"))
                throw std::invalid_argument("Specify exactly one parameters.T or parameters.beta (beta = 0 is infinite temperature)");
            if (p.exists("T") && p["T"].as<double>() <= 0.)
                throw std::invalid_argument("parameters.T must be positive; use beta = 0 for infinite temperature");
            if (!p.exists("THERMALIZATION")) p["THERMALIZATION"] = p["SWEEPS"].as<std::int64_t>() / 10;
            if (!p.exists("S")) p["S"] = p["CONVENTION"].as<std::string>() == "quantum" ? .5 : 1.;
            if (p["MODEL"].as<std::string>() == "Potts" && !p.exists("q"))
                throw std::invalid_argument("Potts requires parameters.q = 3, 4 or 10");
            auto const& e = run.execution;
            if (e.exists("error_variable") != e.exists("error_limit") ||
                (e.exists("error_limit") && e["error_limit"].as<double>() <= 0.))
                throw std::invalid_argument("execution.error_variable and positive error_limit must be supplied together");
            if (e.exists("error_variable")) {
                p["ERROR_VARIABLE"] = e["error_variable"]; p["ERROR_LIMIT"] = e["error_limit"];
            }
            p["PRINT_SWEEPS"] = e["print_sweeps"];
        }, native_mc::spin_output(spinmc::derive));
}
