// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#include "maxent_schema.hpp"
#include <algorithm>
#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/vector.hpp>
#include <alps/maxent.hpp>
#include <alps/run_config.hpp>
#include <cmath>
#include <fstream>
#include <limits>
namespace alps::maxent {
std::string_view schema() { return maxent_schema; }
namespace {
void require(bool ok, const std::string &message) {
    if (!ok)
        throw std::invalid_argument("MaxEnt: " + message);
}
void finite(const std::vector<double> &values, const std::string &name) {
    require(std::all_of(values.begin(), values.end(), [](double x) { return std::isfinite(x); }),
            name + " must be finite");
}
} // namespace
data read_data(const params &supplied) {
    const auto input = resolve_parameters(supplied, schema(), "input");
    data out;
    for (const auto &item : std::vector<std::pair<std::string, std::vector<double> *>>{
             {"values", &out.values},
             {"errors", &out.errors},
             {"covariance", &out.covariance},
             {"tau", &out.tau},
             {"prior_omega", &out.prior_omega},
             {"prior_density", &out.prior_density}})
        if (input.exists(item.first))
            *item.second = input[item.first].as<std::vector<double>>();
    for (const auto *key : {"covariance_dataset", "tau_dataset"})
        if (input.exists(key))
            require(input.exists("data") && input["format"].as<std::string>() == "hdf5",
                    std::string(key) + " requires HDF5 data");
    if (input.exists("data")) {
        require(!input.exists("values") && !input.exists("errors"),
                "choose file data or inline values/errors");
        const auto path = input["data"].as<std::string>();
        if (input["format"].as<std::string>() == "hdf5") {
            hdf5::archive ar(path, "r");
            ar[input["values_dataset"].as<std::string>()] >> out.values;
            if (input.exists("covariance_dataset")) {
                require(!input.exists("covariance"), "covariance supplied twice");
                ar[input["covariance_dataset"].as<std::string>()] >> out.covariance;
            } else if (out.covariance.empty() && !input.exists("covariance_file"))
                ar[input["errors_dataset"].as<std::string>()] >> out.errors;
            if (input.exists("tau_dataset")) {
                require(!input.exists("tau"), "tau supplied twice");
                ar[input["tau_dataset"].as<std::string>()] >> out.tau;
            }
        } else {
            std::ifstream stream(path);
            require(bool(stream), "cannot open data file: " + path);
            int index;
            double value, error;
            while (stream >> std::ws && !stream.eof()) {
                require(bool(stream >> index >> value >> error), "malformed text data: " + path);
                require(index == static_cast<int>(out.values.size()),
                        "text data indices must be consecutive from zero");
                out.values.push_back(value);
                out.errors.push_back(error);
            }
            require(stream.eof(), "malformed text data: " + path);
        }
    }
    if (input.exists("covariance_file")) {
        require(out.covariance.empty(), "covariance supplied twice");
        std::ifstream stream(input["covariance_file"].as<std::string>());
        require(bool(stream), "cannot open covariance file");
        const auto n = out.values.size();
        require(n && n <= std::numeric_limits<std::size_t>::max() / n,
                "invalid covariance dimensions");
        out.covariance.assign(n * n, 0.);
        std::vector<bool> seen(n * n, false);
        int i, j;
        double value;
        while (stream >> std::ws && !stream.eof()) {
            require(bool(stream >> i >> j >> value), "malformed covariance row");
            require(i >= 0 && j >= 0 && static_cast<std::size_t>(i) < n &&
                        static_cast<std::size_t>(j) < n,
                    "covariance index out of range");
            require(!seen[i * n + j], "duplicate covariance entry");
            seen[i * n + j] = true;
            out.covariance[i * n + j] = value;
        }
        require(stream.eof() && std::all_of(seen.begin(), seen.end(), [](bool x) { return x; }),
                "incomplete covariance matrix");
    }
    if (input.exists("prior")) {
        require(out.prior_omega.empty() && out.prior_density.empty(), "prior supplied twice");
        std::ifstream stream(input["prior"].as<std::string>());
        require(bool(stream), "cannot open prior file");
        double omega, density;
        while (stream >> std::ws && !stream.eof()) {
            require(bool(stream >> omega >> density), "malformed prior row");
            out.prior_omega.push_back(omega);
            out.prior_density.push_back(density);
        }
        require(stream.eof(), "malformed prior file");
    }
    return out;
}
params prepare(const params &supplied, const data &input) {
    auto p = resolve_parameters(supplied, schema());
    const auto n = input.values.size();
    require(n >= 4 && n <= static_cast<std::size_t>(std::numeric_limits<int>::max()),
            "requires between 4 and INT_MAX measurements");
    finite(input.values, "data");
    finite(input.errors, "errors");
    finite(input.covariance, "covariance");
    finite(input.tau, "tau");
    require(input.tau.empty() || input.tau.size() == n, "tau must have one value per measurement");
    if (input.covariance.empty()) {
        require(input.errors.size() == n, "errors must have one value per measurement");
        require(
            std::all_of(input.errors.begin(), input.errors.end(), [](double x) { return x > 0; }),
            "errors must be positive");
    } else {
        require(input.covariance.size() == n * n,
                "covariance must be a square matrix matching the data");
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = i + 1; j < n; ++j) {
                const auto a = input.covariance[i * n + j], b = input.covariance[j * n + i];
                require(std::abs(a - b) <= 32 * std::numeric_limits<double>::epsilon() *
                                               std::max({std::abs(a), std::abs(b),
                                                         std::numeric_limits<double>::min()}),
                        "covariance must be symmetric");
            }
    }
    require(p.exists("T") || p.exists("BETA"), "provide T or BETA");
    if (p.exists("T")) {
        require(p["T"].as<double>() > 0, "T must be positive");
        const double beta = 1. / p["T"].as<double>();
        if (p.exists("BETA"))
            require(std::abs(beta - p["BETA"].as<double>()) <= 1e-12 * std::abs(beta),
                    "T and BETA are inconsistent");
        else
            p["BETA"] = beta;
    } else {
        require(p["BETA"].as<double>() > 0, "BETA must be positive");
        p["T"] = 1. / p["BETA"].as<double>();
    }
    for (const auto *key : {"T", "BETA"})
        require(std::isfinite(p[key].as<double>()) && p[key].as<double>() > 0,
                std::string(key) + " must have a finite positive reciprocal");
    if (!p.exists("OMEGA_MIN"))
        p["OMEGA_MIN"] = -p["OMEGA_MAX"].as<double>();
    require(p["OMEGA_MIN"].as<double>() < p["OMEGA_MAX"].as<double>(),
            "OMEGA_MIN must be below OMEGA_MAX");
    for (const auto *key : {"NORM", "ALPHA_MIN", "ALPHA_MAX", "BLOW_UP"})
        require(p[key].as<double>() > 0, std::string(key) + " must be positive");
    require(p["ALPHA_MIN"].as<double>() < p["ALPHA_MAX"].as<double>(),
            "ALPHA_MIN must be below ALPHA_MAX");
    const auto grid = p["FREQUENCY_GRID"].as<std::string>();
    const auto nfreq = p["NFREQ"].as<int>();
    if (grid == "Lorentzian" || grid == "half Lorentzian")
        require(p["CUT"].as<double>() > 0 && p["CUT"].as<double>() < 0.5,
                "CUT must lie between 0 and 0.5");
    if (grid == "quadratic")
        require(p["SPREAD"].as<double>() >= 1, "SPREAD must be at least 1");
    if (grid == "log")
        require(nfreq % 2 == 0 && p["LOG_MIN"].as<double>() > 0 && p["LOG_MIN"].as<double>() < 0.5,
                "log grid needs even NFREQ and 0 < LOG_MIN < 0.5");
    if (p["ENFORCE_NORMALIZATION"].as<bool>())
        require(p.exists("SIGMA_NORMALIZATION") && p["SIGMA_NORMALIZATION"].as<double>() > 0,
                "SIGMA_NORMALIZATION must be positive");
    const auto kernel = p["KERNEL"].as<std::string>(), space = p["DATASPACE"].as<std::string>();
    require(kernel != "anomalous" || space == "frequency",
            "anomalous kernel requires frequency data");
    if (kernel == "Boris")
        require(space == "time" && input.tau.size() == n,
                "Boris kernel needs time data and explicit tau");
    if (space == "frequency" && !p["PARTICLE_HOLE_SYMMETRY"].as<bool>())
        require(n % 2 == 0, "complex frequency data needs paired real/imaginary values");
    const auto prior = p["DEFAULT_MODEL"].as<std::string>();
    auto positive = [&](const char *key) {
        require(p.exists(key) && p[key].as<double>() > 0,
                std::string(key) + " must be supplied and positive");
    };
    if (prior.find("gaussian") != std::string::npos && prior != "twogaussians")
        positive("SIGMA");
    if (prior == "twogaussians") {
        positive("SIGMA1");
        positive("SIGMA2");
        require(p.exists("SHIFT2"), "twogaussians requires SHIFT2");
        require(p["NORM1"].as<double>() >= 0 && p["NORM1"].as<double>() <= 1,
                "NORM1 must lie in [0,1]");
    }
    if (prior == "shifted gaussian" || prior == "double gaussian" ||
        prior == "general double gaussian")
        require(p.exists("SHIFT"), "shifted Gaussian prior requires SHIFT");
    if (prior == "general double gaussian")
        positive("BOSE_NORM");
    if (prior == "linear rise exp decay" || prior == "quadratic rise exp decay") {
        positive("LAMBDA");
        require(p["OMEGA_MIN"].as<double>() >= 0,
                "exponential prior requires nonnegative frequencies");
    }
    if (prior == "tabulated") {
        finite(input.prior_omega, "prior coordinates");
        finite(input.prior_density, "prior density");
        require(input.prior_omega.size() >= 2 &&
                    input.prior_omega.size() == input.prior_density.size(),
                "tabulated prior needs matching coordinates and density");
        require(input.prior_omega.front() <= p["OMEGA_MIN"].as<double>() &&
                    input.prior_omega.back() >= p["OMEGA_MAX"].as<double>(),
                "tabulated prior must cover the frequency range");
        require(std::any_of(input.prior_density.begin(), input.prior_density.end(),
                            [](double x) { return x > 0; }),
                "prior density cannot be identically zero");
        for (std::size_t i = 0; i < input.prior_omega.size(); ++i) {
            require(input.prior_density[i] >= 0, "prior density must be nonnegative");
            if (i)
                require(input.prior_omega[i] > input.prior_omega[i - 1],
                        "prior coordinates must increase");
        }
    }
    return p;
}
data prepare_run(run_configuration &run) {
    run = resolve_run_configuration(run, schema(), std::filesystem::current_path());
    auto input = read_data(run.input);
    run.parameters = prepare(run.parameters, input);
    for (const auto *key : {"T", "BETA", "OMEGA_MIN"})
        if (!run.origins.count(std::string("parameters.") + key))
            run.origins[std::string("parameters.") + key] = "derived";
    return input;
}
} // namespace alps::maxent

bool alps::solvers::maxent(const alps::run_configuration &supplied) {
    auto run = supplied;
    auto data = alps::maxent::prepare_run(run);
    const auto output = run.output["results"].as<std::string>();
    const bool completed = maxent(run.parameters, data, output,
                                 run.execution["time_limit"].as<int>(),
                                 run.output["text"].as<bool>());
    if (completed) {
        alps::hdf5::archive ar(output, "a");
        ar["/run_config"] << run;
    }
    return completed;
}
