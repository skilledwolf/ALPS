// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
#include <alps/params.hpp>
#include <alps/maxent.hpp>

#include <stdexcept>
#include <string>

// Calling the public entry point pulls the complete static solver into this
// executable. Link only ALPS::maxent to verify its exported dependency closure.
int main() {
  alps::params parameters;
  parameters["OMEGA_MAX"] = 4.0;
  parameters["BETA"] = 2.0;
  parameters["NFREQ"] = 20;
  parameters["DEFAULT_MODEL"] = std::string("flat");
  try {
    alps::solvers::maxent(parameters, {}, "unused.out.h5");
  } catch (std::invalid_argument const& error) {
    return std::string(error.what()) .find("measurements") != std::string::npos ? 0 : 1;
  }
  return 2;
}
