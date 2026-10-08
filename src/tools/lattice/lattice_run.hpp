// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#pragma once
#include <alps/parameter.h>
#include <alps/parser/xslt_path.h>
#include <alps/run_config.hpp>
#include <filesystem>
#include <stdexcept>
#include <string>

// The lattice tools read the [parameters] of a TOML run file: exactly one
// LATTICE or GRAPH and its lattice parameters. As for the lattice-model
// applications, the lattice library is input.lattice_library.
inline alps::Parameters lattice_run(std::filesystem::path const& file, std::string const& application) {
    const std::string base = "application = \"" + application + "\"\nschema_version = 1\n"
        "[parameters.LATTICE]\ntype = \"string\"\n[parameters.GRAPH]\ntype = \"string\"\n"
        "[input.lattice_library]\ntype = \"path\"\n[output]\n[execution]\n";
    const auto run = alps::load_run_configuration(file, alps::extend_run_schema(file, base,
        [](std::string const& key) -> char const* {
            if (key == "LATTICE_LIBRARY")
                throw std::invalid_argument("Retired parameter LATTICE_LIBRARY; use input.lattice_library");
            return nullptr;
        }));
    if (run.parameters.exists("LATTICE") == run.parameters.exists("GRAPH"))
        throw std::invalid_argument("Specify exactly one parameters.LATTICE or parameters.GRAPH");
    alps::Parameters parameters(run.parameters);
    parameters["LATTICE_LIBRARY"] =
        alps::search_xml_library_path(run.input.value_or<std::string>("lattice_library", "lattices.xml"));
    return parameters;
}
