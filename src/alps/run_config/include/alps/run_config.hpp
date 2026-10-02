// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#pragma once
#include <alps/params.hpp>
#include <alps/run_config_export.h>
#include <filesystem>
#include <map>
#include <string_view>
namespace alps {
// Application schemas are TOML supplied by the application. Parser types stay
// private; native and language-binding callers share the same validation.
struct ALPS_RUN_CONFIG_DECL run_configuration {
    std::string application;
    int schema_version = 1;
    params parameters, input, output, execution;
    std::map<std::string, std::string> origins;
    void save(hdf5::archive &) const;
};
ALPS_RUN_CONFIG_DECL run_configuration load_run_configuration(const std::filesystem::path &filename,
                                                              std::string_view schema);
ALPS_RUN_CONFIG_DECL params resolve_parameters(const params &supplied, std::string_view schema,
                                               const std::string &section = "parameters");
} // namespace alps
