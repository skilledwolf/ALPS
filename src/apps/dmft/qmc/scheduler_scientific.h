// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#pragma once
#include <alps/parameter.h>
#include <alps/params.hpp>
// The multiband scheduler CT-INT still receives Parameters from MCRun. Its
// numerical helpers use typed params. This private snapshot disappears with that
// scheduler constructor; it is not a runtime input parser.
inline alps::params scheduler_scientific(const alps::Parameters& old) {
  alps::params result;
  for (const char* key : {"BETA", "MU", "H", "U", "J", "U'"})
    if (old.defined(key)) result[key] = static_cast<double>(old[key]);
  result["FLAVORS"] = static_cast<int>(old.value_or_default("FLAVORS", 2));
  result["SITES"] = static_cast<int>(old.value_or_default("SITES", 1));
  for (int f = 0; f < result["FLAVORS"].as<int>(); ++f)
    for (const char* prefix : {"EPS_", "EPSSQ_"}) {
      const auto key = std::string(prefix) + std::to_string(f);
      if (old.defined(key)) result[key] = static_cast<double>(old[key]);
    }
  return result;
}
inline alps::params scheduler_input(const alps::Parameters& old) {
  alps::params input;
  if (old.defined("U_MATRIX")) input["interaction_matrix"] = static_cast<std::string>(old["U_MATRIX"]);
  return input;
}
