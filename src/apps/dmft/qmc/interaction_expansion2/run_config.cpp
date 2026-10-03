// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#include <alps/ctint.hpp>
#include "run_config.hpp"
#include "ctint_schema.hpp"
#include <alps/hdf5/complex.hpp>
#include <alps/hdf5/vector.hpp>
#include <boost/math/constants/constants.hpp>
#include <cmath>
#include <filesystem>
#include <limits>

std::string_view alps::ctint::schema() { return ctint_schema; }

alps::params alps::ctint::prepare_parameters(const params &supplied) {
    auto parameters = resolve_parameters(supplied, schema());
    if (parameters["BETA"].as<double>() <= 0.0)
        throw std::invalid_argument("CT-INT BETA must be positive");
    if (!parameters.exists("NMATSUBARA_MEASUREMENTS"))
        parameters["NMATSUBARA_MEASUREMENTS"] = parameters["NMATSUBARA"];
    if (!parameters.exists("NSELF"))
        parameters["NSELF"] = 10 * parameters["N"].as<std::int64_t>();
    if (parameters["NMATSUBARA_MEASUREMENTS"].as<int>() > parameters["NMATSUBARA"].as<int>())
        throw std::invalid_argument("CT-INT NMATSUBARA_MEASUREMENTS must not exceed NMATSUBARA");
    const auto sweeps = parameters["SWEEPS"].as<std::uint64_t>();
    const auto thermalization = parameters["THERMALIZATION"].as<std::uint64_t>();
    const auto period = parameters["MEASUREMENT_PERIOD"].as<std::uint64_t>();
    if (sweeps > std::numeric_limits<std::uint64_t>::max() - thermalization ||
        sweeps + thermalization > std::numeric_limits<std::uint64_t>::max() - period)
        throw std::invalid_argument("CT-INT total sweep count exceeds the counter range");
    if (parameters.exists("TAU_DISCRETIZATION_FOR_EXP")) {
        const auto frequencies = parameters["NMATSUBARA"].as<std::size_t>();
        const auto slices = parameters["TAU_DISCRETIZATION_FOR_EXP"].as<std::size_t>();
        if (frequencies > std::numeric_limits<std::size_t>::max() /
                          (2 * sizeof(std::complex<double>)) / slices)
            throw std::invalid_argument("CT-INT exponential table dimensions overflow");
    }
    return parameters;
}

void read_ctint_bare_green(const alps::params &parameters, const alps::params &input,
                          matsubara_green_function_t &green) {
    const auto n = parameters["NMATSUBARA"].as<std::size_t>();
    if (input["atomic"].as<bool>()) {
        const auto beta = parameters["BETA"].as<double>();
        for (unsigned int flavor = 0; flavor < green.nflavor(); ++flavor)
            for (std::size_t frequency = 0; frequency < n; ++frequency)
                green(frequency, 0, 0, flavor) = {0.0, -beta / ((2.0 * frequency + 1.0) *
                    boost::math::constants::pi<double>())};
        return;
    }
    alps::hdf5::archive archive(input["g0"].as<std::string>(), "r");
    read_flavor_vectors(archive, "/G0", green);
}

void alps::ctint::prepare_run(run_configuration &run) {
    run = resolve_run_configuration(run, schema(), std::filesystem::current_path());
    run.parameters = prepare_parameters(run.parameters);
    for (const auto *key : {"NMATSUBARA_MEASUREMENTS", "NSELF"})
        if (!run.origins.count(std::string("parameters.") + key))
            run.origins[std::string("parameters.") + key] = "derived";
    if (run.execution["check_interval"].as<double>() <= 0.0)
        throw std::invalid_argument("CT-INT execution.check_interval must be positive");
    const bool atomic = run.input["atomic"].as<bool>();
    if (atomic == run.input.exists("g0"))
        throw std::invalid_argument("CT-INT requires either input.g0 or input.atomic=true");
    if (atomic && (run.parameters["MU"].as<double>() != 0.0 ||
                   run.parameters["H"].as<double>() != 0.0))
        throw std::invalid_argument("CT-INT atomic input requires MU=0 and H=0");
    matsubara_green_function_t green(run.parameters["NMATSUBARA"].as<unsigned int>(), 1, 2);
    read_ctint_bare_green(run.parameters, run.input, green);
}
