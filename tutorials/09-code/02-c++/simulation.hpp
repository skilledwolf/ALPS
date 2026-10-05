/*****************************************************************************
 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations
 *
 * ALPS Libraries
 *
 * Copyright (C) 2003 by Brigitte Surer
 *                       and Jan Gukelberger
 *
 * ALPS Project: https://alps.comp-phys.org/
 * SPDX-License-Identifier: MIT
 *
 *****************************************************************************/

#pragma once
#include <alps/alea/batch.hpp>
#include <alps/alea/autocorr.hpp>
#include <alps/alea/hdf5.hpp>
#include <alps/alea/transform.hpp>
#include <alps/alea/transformer.hpp>
#include <alps/run_config.hpp>
#include <boost/random.hpp>
#include <array>
#include <cmath>
#include <iostream>
#include <set>
#include <sstream>
#include <vector>

inline constexpr char schema[] = R"toml(application = "cpp-ising"
schema_version = 1
[parameters.L]
type = "int64"
default = 16
min = 2
max = 2147483647
[parameters.BETA]
type = "float64"
default = 0.4
min = 0.0
[parameters.THERMALIZATION]
type = "int64"
default = 2500
min = 0
[parameters.SWEEPS]
type = "int64"
default = 5000
min = 2
[input]
[output.results]
type = "path"
required = true
[execution.seed]
type = "int64"
default = 42
min = 0
max = 2147483647
[execution.bins]
type = "int64"
default = 64
min = 2
)toml";

class Simulation {
public:
    explicit Simulation(alps::run_configuration const& run)
        : L(run.parameters["L"].as<int>()), beta(run.parameters["BETA"].as<double>()),
          random(boost::mt19937(run.execution["seed"].as<int>()), boost::uniform_real<>()),
          spins(L, std::vector<int>(L)), samples(5, run.execution["bins"].as<size_t>()) {
        for (auto& row : spins) for (auto& spin : row) spin = 2*randint(2)-1;
    }
    void step();
    std::vector<double> observables() const;
    void run(alps::run_configuration const& config) {
        for (int64_t i=0; i<config.parameters["THERMALIZATION"].as<int64_t>(); ++i) step();
        for (int64_t i=0; i<config.parameters["SWEEPS"].as<int64_t>(); ++i) {
            step();
            auto values = observables();
            samples << alps::alea::make_adapter(values);
            for (size_t j=0; j<names.size(); ++j) diagnostics[j] << alps::alea::make_adapter(values[j]);
        }
        const auto joint = samples.result();
        alps::hdf5::save_checkpoint(config.output["results"].as<std::string>(), [&](auto& archive) {
            archive["/parameters"] << config.parameters;
            archive["/run_config"] << config;
            alps::alea::hdf5_serializer raw(archive, "/simulation");
            serialize(raw, "joint", joint);
            alps::alea::hdf5_serializer results(archive, "/simulation/results");
            alps::alea::hdf5_serializer analysis(archive, "/simulation/realizations/0/clones/0/autocorrelation");
            for (size_t j=0; j<names.size(); ++j) {
                Eigen::Matrix<double,1,5> select = Eigen::Matrix<double,1,5>::Zero();
                select(j) = 1;
                auto result = alps::alea::transform(alps::alea::jackknife_prop(),
                    alps::alea::linear_transformer<double>(select), joint);
                const auto diagnostic = diagnostics[j].result();
                const auto name = archive.encode_segment(names[j]);
                serialize(results, name, result);
                serialize(analysis, name, diagnostic);
                std::cout << names[j] << ": " << result.mean()(0) << " +/- " << result.stderror()(0)
                          << "; convergence (0=yes, 1=unknown, 2=no): " << diagnostic.converged_errors()(0);
                if (diagnostic.tau_available()) std::cout << "; tau = " << diagnostic.tau()(0);
                std::cout << '\n';
            }
        });
    }
private:
    int randint(int maximum) { return int(maximum*random()); }
    int L;
    double beta;
    boost::variate_generator<boost::mt19937, boost::uniform_real<>> random;
    std::vector<std::vector<int>> spins;
    alps::alea::batch_acc<double> samples;
    std::array<alps::alea::autocorr_acc<double>,5> diagnostics;
    inline static const std::array<std::string,5> names{"E", "m", "|m|", "m^2", "m^4"};
};

inline int run_tutorial(int argc, char** argv) {
    try {
        bool validate = false;
        std::vector<alps::run_configuration> runs;
        std::set<std::filesystem::path> inputs, outputs;
        for (int i=1; i<argc; ++i) {
            const std::string arg(argv[i]);
            if (arg == "--schema") { std::cout << schema; return 0; }
            if (arg == "--help" || arg == "-h") {
                std::cout << "Usage: ising [--validate] run.toml [run.toml ...] | ising --schema\n"
                          << "No arguments runs the tutorial beta scan.\n";
                return 0;
            }
            if (arg == "--validate") { validate = true; continue; }
            if (arg.empty() || arg.front() == '-') throw std::invalid_argument("Unknown option: " + arg);
            inputs.insert(std::filesystem::weakly_canonical(arg));
            runs.push_back(alps::load_run_configuration(arg, schema));
        }
        if (validate && runs.empty()) throw std::invalid_argument("--validate requires a run file");
        if (runs.empty()) for (int i=0; i<=10; ++i) {
            alps::run_configuration run;
            run.parameters["BETA"] = i/10.;
            std::ostringstream filename;
            filename << "ising.L_16beta_" << i/10. << ".h5";
            run.output["results"] = filename.str();
            runs.push_back(alps::resolve_run_configuration(run, schema));
        }
        for (auto const& run : runs) {
            if (run.execution["bins"].as<int64_t>() % 2)
                throw std::invalid_argument("execution.bins must be even");
            auto output = std::filesystem::weakly_canonical(run.output["results"].as<std::string>());
            if (inputs.count(output) || !outputs.insert(output).second)
                throw std::invalid_argument("Outputs must be distinct and must not replace run files");
        }
        for (auto const& run : runs) {
            if (validate) std::cout << "Valid Ising configuration\n";
            else Simulation(run).run(run);
        }
        return 0;
    } catch (std::exception const& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
