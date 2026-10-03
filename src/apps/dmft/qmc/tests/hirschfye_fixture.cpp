// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#include <alps/run_config.hpp>
#include <alps/hdf5/complex.hpp>
#include <alps/hdf5/vector.hpp>
#include <cmath>
#include <complex>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr unsigned int frequencies = 4;
constexpr unsigned int time_points = 5;
constexpr double beta = 2.0;

std::complex<double> free_green(unsigned int frequency) {
  return {0.0, -beta / ((2.0 * frequency + 1.0) * std::acos(-1.0))};
}

void verify_output(const std::string& filename) {
  alps::hdf5::archive archive(filename, "r");
  alps::run_configuration run;
  archive["/run_config"] >> run;
  if (run.application != "hirschfye" || run.schema_version != 1 ||
      run.parameters["BETA"].as<double>() != beta ||
      run.parameters["U"].as<double>() != 0.0 ||
      run.parameters["MU"].as<double>() != 0.0 ||
      run.parameters["EPSSQ_0"].as<double>() != 0.0 ||
      run.parameters["EPSSQ_1"].as<double>() != 0.0 ||
      run.parameters["N"].as<unsigned int>() != time_points - 1 ||
      run.parameters["NMATSUBARA"].as<unsigned int>() != frequencies ||
      !run.input.exists("g0") || !run.output.exists("results"))
    throw std::runtime_error("Hirsch-Fye result does not retain the resolved typed run metadata");
  for (unsigned int flavor = 0; flavor < 2; ++flavor) {
    std::vector<double> tau;
    std::vector<std::complex<double>> omega;
    const auto channel = "/" + std::to_string(flavor) + "/mean/value";
    archive["/G_tau" + channel] >> tau;
    archive["/G_omega" + channel] >> omega;
    if (tau.size() != time_points || omega.size() != frequencies)
      throw std::runtime_error("Hirsch-Fye free result has incorrect dimensions");
    for (const auto value : tau)
      if (!std::isfinite(value) || std::abs(value + 0.5) > 1.e-10)
        throw std::runtime_error("Hirsch-Fye U=0 G(tau) differs from the exact -1/2 result");
    for (unsigned int frequency = 0; frequency < frequencies; ++frequency)
      if (!std::isfinite(omega[frequency].real()) || !std::isfinite(omega[frequency].imag()) ||
          std::abs(omega[frequency] - free_green(frequency)) > 1.e-10)
        throw std::runtime_error("Hirsch-Fye U=0 G(iw) differs from the exact 1/iw result");
  }
}

void write_fixture(const std::string& filename, const std::string& kind) {
  if (kind != "valid" && kind != "short" && kind != "oversized" && kind != "real" &&
      kind != "missing-flavor" && kind != "nonfinite")
    throw std::invalid_argument("Unknown Hirsch-Fye fixture kind: " + kind);
  alps::hdf5::archive archive(filename, "w");
  for (unsigned int flavor = 0; flavor < (kind == "missing-flavor" ? 1u : 2u); ++flavor) {
    const auto path = "/G0_" + std::to_string(flavor);
    std::vector<std::complex<double>> values(
        kind == "oversized" ? frequencies + 1 : kind == "short" ? frequencies - 1 : frequencies);
    for (unsigned int frequency = 0; frequency < values.size(); ++frequency)
      values[frequency] = free_green(frequency);
    if (kind == "nonfinite") values[0] = {0.0, std::numeric_limits<double>::infinity()};
    if (kind == "real") archive[path] << std::vector<double>(frequencies, -0.5);
    else archive[path] << values;
  }
}
}

int main(int argc, char** argv) {
  try {
    if (argc != 3) throw std::invalid_argument("Usage: hirschfye_fixture output.h5 kind");
    if (std::string(argv[2]) == "verify-output") verify_output(argv[1]);
    else write_fixture(argv[1], argv[2]);
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
