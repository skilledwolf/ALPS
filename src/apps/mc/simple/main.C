// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#include "simulation.h"
#include "schema.hpp"
#include "../native_driver.hpp"

int main(int argc, char** argv) {
    return native_mc::main<simplemc::simulation>(argc, argv, "simplemc", simplemc_schema,
        [](std::string const& key, toml::node const&) -> char const* {
            return key.size() > 1 && key[0] == 'J' &&
                key.find_first_not_of("0123456789", 1) == std::string::npos ? "float64" : nullptr;
        },
        [](alps::params& p, alps::run_configuration const&) {
            if (!p.exists("THERMALIZATION")) p["THERMALIZATION"] = p["SWEEPS"].as<std::int64_t>() / 8;
            if (p.exists("T") && p["T"].as<double>() <= 0.)
                throw std::invalid_argument("parameters.T must be positive; omit it for infinite temperature");
        }, simplemc::simulation::derive,
        [](simplemc::simulation const& simulation, std::filesystem::path const& path) { simulation.snapshot(path); });
}
