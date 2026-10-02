// Copyright (C) 2010 Sebastian Fuchs, Thomas Pruschke, Matthias Troyer.
// Modifications (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#include <alps/hdf5/archive.hpp>
#include <alps/maxent.hpp>
#include <alps/run_config.hpp>
#include <iostream>
#include <string>
int main(int argc, char **argv) {
    try {
        bool validate = false;
        std::string file;
        for (int i = 1; i < argc; ++i) {
            const std::string arg = argv[i];
            if (arg == "--help" || arg == "-h") {
                std::cout << "Usage: maxent [--validate] run.toml\n"
                          << "Input, output and execution settings belong in the TOML run file.\n";
                return 0;
            }
        }
        for (int i = 1; i < argc; ++i) {
            const std::string arg = argv[i];
            if (arg == "--validate")
                validate = true;
            else if (arg.empty() || arg[0] == '-')
                throw std::invalid_argument("unknown option: " + arg);
            else if (file.empty())
                file = arg;
            else
                throw std::invalid_argument("expected one TOML run file");
        }
        if (file.empty())
            throw std::invalid_argument("No TOML run file specified");
        auto run = alps::load_run_configuration(file, alps::maxent::schema());
        const auto data = alps::maxent::read_data(run.input);
        run.parameters = alps::maxent::prepare(run.parameters, data);
        for (const auto *key : {"T", "BETA", "OMEGA_MIN"})
            if (!run.origins.count(std::string("parameters.") + key))
                run.origins[std::string("parameters.") + key] = "derived";
        const auto output = run.output["results"].as<std::string>();
        for (const auto *key : {"data", "covariance_file", "prior"})
            if (run.input.exists(key)) {
                const auto source = std::filesystem::path(run.input[key].as<std::string>());
                if (std::filesystem::weakly_canonical(source) ==
                    std::filesystem::weakly_canonical(output))
                    throw std::invalid_argument("output must not replace an input file");
            }
        if (std::filesystem::weakly_canonical(file) == std::filesystem::weakly_canonical(output))
            throw std::invalid_argument("output must not replace the TOML run file");
        if (validate) {
            std::cout << "Valid MaxEnt configuration: " << file << '\n';
            return 0;
        }
        const bool completed = alps::solvers::maxent(run.parameters, data, output,
                                                     run.execution["time_limit"].as<int>(),
                                                     run.output["text"].as<bool>());
        // A stopped run without a result has no configuration archive to annotate.
        if (completed) {
            alps::hdf5::archive ar(output, "a");
            ar["/run_config"] << run;
        }
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "maxent: " << error.what() << '\n';
        return 1;
    }
}
