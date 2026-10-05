// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#pragma once
#include <alps/params.hpp>
#include <alps/run_config_export.h>
#include <filesystem>
#include <functional>
#include <map>
#include <string_view>
namespace alps {
// Application schemas are TOML supplied by the application. Parser types stay
// private; native and language-binding callers share the same validation.
// Identity and version come from the schema, not the user-written run file.
struct ALPS_RUN_CONFIG_DECL run_configuration {
    std::string application;
    int schema_version = 1;
    params parameters, input, output, execution;
    std::map<std::string, std::string> origins;
    // Transient source path for diagnostics and protection against overwriting
    // the loaded run file. Empty for programmatic or archive-loaded runs.
    std::string source_file;
    void save(hdf5::archive &) const;
    void load(hdf5::archive &);
};
// Add undeclared model parameters found in a run file. The callback may select
// a schema type for a key (e.g. a coupling array); nullptr infers scalar types.
// An empty path only normalizes the base schema. TOML parser types stay private.
ALPS_RUN_CONFIG_DECL std::string extend_run_schema(const std::filesystem::path &file,
    std::string_view base, const std::function<char const*(std::string const&)> &type = {});
// Programmatic callers use the same section rules and defaults as TOML files.
// Relative paths are resolved only when a base directory is supplied. Loading
// and resolution reject path outputs that would replace the run file, an input
// or another output.
ALPS_RUN_CONFIG_DECL run_configuration resolve_run_configuration(
    const run_configuration &supplied, std::string_view schema,
    const std::filesystem::path &base_directory = {});
ALPS_RUN_CONFIG_DECL run_configuration load_run_configuration(const std::filesystem::path &filename,
                                                              std::string_view schema);
// Serialize the four sections without schema metadata. Values must be finite
// and representable in TOML; schema validation remains application-owned.
ALPS_RUN_CONFIG_DECL std::string format_run_configuration(const run_configuration &run);
ALPS_RUN_CONFIG_DECL params resolve_parameters(const params &supplied, std::string_view schema,
                                               const std::string &section = "parameters");
// Project a parent run onto a child application's declared keys. Resolution
// subsequently validates values and inserts the child's defaults.
ALPS_RUN_CONFIG_DECL params select_parameters(const params &supplied, std::string_view schema,
                                              const std::string &section = "parameters");
} // namespace alps
