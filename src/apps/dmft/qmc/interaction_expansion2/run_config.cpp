// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#include <alps/ctint.hpp>
#include "run_config.hpp"
#include "ctint_schema.hpp"
#include "U_matrix.h"
#include <alps/hdf5/complex.hpp>
#include <alps/hdf5/vector.hpp>
#include <boost/math/constants/constants.hpp>
#include <cmath>
#include <filesystem>
#include <limits>
#include <sstream>
#include <toml++/toml.hpp>

std::string alps::ctint::schema(const params &parameters) {
    const auto flavors = parameters.value_or<std::int64_t>("FLAVORS", 2);
    if (flavors < 2 || flavors > 128)
        throw std::invalid_argument("CT-INT FLAVORS must be in [2,128]");
    std::ostringstream out;
    out << ctint_schema;
    for (int flavor = 2; flavor < flavors; ++flavor)
        for (const auto *prefix : {"EPS_", "EPSSQ_"})
            out << "\n[parameters." << prefix << flavor << "]\ntype = \"float64\"\ndefault = "
                << (prefix == std::string_view("EPS_") ? "0.0" : "1.0") << '\n';
    return out.str();
}

std::string alps::ctint::schema_for_run(const std::filesystem::path &file) {
    if (file.empty()) return schema();
    const auto raw = toml::parse_file(file.string());
    params selectors;
    if (const auto *node = raw["parameters"]["FLAVORS"].node()) {
        if (const auto value = node->value_exact<std::int64_t>()) selectors["FLAVORS"] = *value;
        else throw std::invalid_argument("CT-INT FLAVORS must be an integer");
    }
    return schema(selectors);
}

alps::params alps::ctint::prepare_parameters(const params &supplied) {
    auto parameters = resolve_parameters(supplied, schema(supplied));
    const auto beta = parameters["BETA"].as<double>();
    if (beta <= 0.0)
        throw std::invalid_argument("CT-INT BETA must be positive");
    if (!parameters.exists("NMATSUBARA_MEASUREMENTS"))
        parameters["NMATSUBARA_MEASUREMENTS"] = parameters["NMATSUBARA"];
    if (!parameters.exists("NSELF"))
        parameters["NSELF"] = 10 * parameters["N"].as<std::int64_t>();
    if (parameters["NMATSUBARA_MEASUREMENTS"].as<int>() > parameters["NMATSUBARA"].as<int>())
        throw std::invalid_argument("CT-INT NMATSUBARA_MEASUREMENTS must not exceed NMATSUBARA");
    const auto flavors = parameters["FLAVORS"].as<std::uint64_t>();
    if ((parameters["N"].as<std::uint64_t>() + 1) * flavors > std::numeric_limits<unsigned>::max() ||
        parameters["NMATSUBARA"].as<std::uint64_t>() * flavors > std::numeric_limits<unsigned>::max())
        throw std::invalid_argument("CT-INT Green-function dimensions exceed the supported storage range");
    if (!std::isfinite(beta * beta) ||
        !std::isfinite((2. * parameters["NMATSUBARA"].as<double>() - 1.) *
                       boost::math::constants::pi<double>() / beta))
        throw std::invalid_argument("CT-INT BETA and NMATSUBARA exceed the Fourier arithmetic range");
    for (std::uint64_t flavor = 0; flavor < flavors; ++flavor) {
        const auto mu = parameters["MU"].as<double>() + (flavor % 2 ? 1. : -1.) * parameters["H"].as<double>();
        const auto eps = parameters["EPS_" + std::to_string(flavor)].as<double>();
        const auto epssq = parameters["EPSSQ_" + std::to_string(flavor)].as<double>();
        const auto second = eps - mu;
        const auto third = epssq - 2. * mu * eps + mu * mu;
        if (!std::isfinite(mu) || !std::isfinite(second) || !std::isfinite(third) ||
            !std::isfinite(beta * second) || !std::isfinite(beta * beta * third))
            throw std::invalid_argument("CT-INT flavor moments exceed the Fourier arithmetic range");
    }
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

void alps::ctint::validate_interaction(const params &parameters, const params &input) {
    if (!input.exists("interaction_matrix")) {
        if (!parameters.exists("U"))
            throw std::invalid_argument("CT-INT requires parameters.U or input.interaction_matrix");
        if (parameters["FLAVORS"].as<unsigned>() % 2)
            throw std::invalid_argument("CT-INT generated Hund interactions require paired FLAVORS or input.interaction_matrix");
    }
    const U_matrix interaction(parameters, input);
    for (spin_t flavor = 0; flavor < interaction.nf(); ++flavor) {
        double sum = 0., magnitude = 0.;
        if (interaction(flavor, flavor) != 0.)
            throw std::invalid_argument("CT-INT interaction matrix must have a zero diagonal");
        for (spin_t other = 0; other < interaction.nf(); ++other) {
            const auto value = interaction(flavor, other);
            if (value != interaction(other, flavor))
                throw std::invalid_argument("CT-INT interaction matrix must be symmetric");
            sum += value;
            magnitude += std::abs(value);
            if (!std::isfinite(parameters["BETA"].as<double>() * std::abs(value) *
                               (interaction.n_nonzero() / 2.)))
                throw std::invalid_argument("CT-INT interaction weights exceed the arithmetic range");
        }
        const auto mu = parameters["MU"].as<double>() + (flavor % 2 ? 1. : -1.) * parameters["H"].as<double>();
        if (!std::isfinite(sum) || !std::isfinite(magnitude * magnitude) || !std::isfinite(mu + .5 * sum))
            throw std::invalid_argument("CT-INT interaction moments exceed the arithmetic range");
    }
}

void alps::ctint::validate_execution(const params &execution) {
    if (execution["check_interval"].as<double>() <= 0.)
        throw std::invalid_argument("CT-INT execution.check_interval must be positive");
    const auto bins = execution["bins"].as<std::uint64_t>();
    if (bins < 2 || bins % 2)
        throw std::invalid_argument("CT-INT execution.bins must be even and at least two");
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
    run = resolve_run_configuration(run, schema(run.parameters), std::filesystem::current_path());
    run.parameters = prepare_parameters(run.parameters);
    for (const auto *key : {"NMATSUBARA_MEASUREMENTS", "NSELF"})
        if (!run.origins.count(std::string("parameters.") + key))
            run.origins[std::string("parameters.") + key] = "derived";
    validate_execution(run.execution);
    const auto output = std::filesystem::path(run.output["results"].as<std::string>());
    if (std::filesystem::is_directory(output) || !std::filesystem::is_directory(output.parent_path()))
        throw std::invalid_argument("CT-INT output.results must name a file in an existing directory");
    for (const auto *key : {"g0", "interaction_matrix"})
        if (run.input.exists(key) && !std::filesystem::is_regular_file(run.input[key].as<std::string>()))
            throw std::invalid_argument(std::string("CT-INT input.") + key + " must be an existing file");
    validate_interaction(run.parameters, run.input);
    const bool atomic = run.input["atomic"].as<bool>();
    if (atomic == run.input.exists("g0"))
        throw std::invalid_argument("CT-INT requires either input.g0 or input.atomic=true");
    if (atomic && (run.parameters["MU"].as<double>() != 0.0 ||
                   run.parameters["H"].as<double>() != 0.0))
        throw std::invalid_argument("CT-INT atomic input requires MU=0 and H=0");
    if (atomic)
        for (unsigned flavor = 0; flavor < run.parameters["FLAVORS"].as<unsigned>(); ++flavor)
          for (const auto *prefix : {"EPS_", "EPSSQ_"}) {
            const auto key = std::string(prefix) + std::to_string(flavor);
            auto &origin = run.origins.at(std::string("parameters.") + key);
            if (origin == "default" || origin == "derived") {
                run.parameters[key] = 0.0;
                origin = "derived";
            } else if (run.parameters[key].as<double>() != 0.0)
                throw std::invalid_argument("CT-INT atomic input requires " + key + "=0");
        }
    matsubara_green_function_t green(run.parameters["NMATSUBARA"].as<unsigned int>(), 1,
                                     run.parameters["FLAVORS"].as<unsigned>());
    read_ctint_bare_green(run.parameters, run.input, green);
}
