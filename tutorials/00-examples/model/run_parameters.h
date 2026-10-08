// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#pragma once
#include <alps/expression/symbol_table.h>
#include <alps/run_config.hpp>
#include <stdexcept>
#include <string>

// The examples read the model and lattice from the [parameters] of a TOML run
// file, such as parm_numeric.toml.
inline alps::SymbolTable run_parameters(int argc, char** argv) {
  if (argc != 2)
    throw std::invalid_argument(std::string("Usage: ") + argv[0] + " run.toml");
  const char schema[] = "application = \"model-example\"\nschema_version = 1\n"
                        "[input]\n[output]\n[execution]\n";
  return alps::SymbolTable(
    alps::load_run_configuration(argv[1], alps::extend_run_schema(argv[1], schema)).parameters);
}
