// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#pragma once
#include <alps/hdf5/archive.hpp>
#include <alps/parameter.h>
#include <alps/parser/xslt_path.h>
#include <alps/run_config.hpp>
#include <filesystem>
#include <iostream>
#include <set>
#include <string>
#include <vector>

// Applications on the lattice and model libraries read one TOML run per
// calculation and publish its results to output.results.
//
// An application type App provides
//   static void print_copyright(std::ostream&);
//   template <class F> static void with_task(alps::run_configuration const&, F const&);
// where with_task constructs the run's task, which checks its parameters, and
// passes it to F. The task provides run() and save(alps::hdf5::archive&).
namespace lattice_model {
inline std::string schema(std::filesystem::path const& file, char const* base) {
    return alps::extend_run_schema(file, base, [](std::string const& key) -> char const* {
        if (key == "LATTICE_LIBRARY" || key == "MODEL_LIBRARY")
            throw std::invalid_argument("Retired parameter " + key + "; use input.lattice_library or input.model_library");
        return nullptr;
    });
}

inline alps::Parameters parameters(alps::run_configuration const& run) {
    alps::Parameters p(run.parameters);
    p["LATTICE_LIBRARY"] = run.input["lattice_library"].as<std::string>();
    p["MODEL_LIBRARY"] = run.input["model_library"].as<std::string>();
    return p;
}

template <class App>
int main(int argc, char** argv, char const* application, char const* base_schema) {
    try {
        bool validate = false, show_schema = false;
        std::vector<std::filesystem::path> files;
        for (int i = 1; i < argc; ++i) {
            const std::string argument(argv[i]);
            if (argument == "--help" || argument == "-h") {
                std::cout << "Usage: " << application << " [--validate] run.toml [run.toml ...] | "
                          << application << " --schema [run.toml]\n"
                             "Model and lattice parameters belong in [parameters], libraries in\n"
                             "input.lattice_library/model_library and the HDF5 file in output.results.\n";
                return 0;
            }
            if (argument == "--validate") validate = true;
            else if (argument == "--schema") show_schema = true;
            else if (argument.empty() || argument.front() == '-') throw std::invalid_argument("Unknown option: " + argument);
            else files.emplace_back(argument);
        }
        if (show_schema) {
            if (files.size() > 1 || validate) throw std::invalid_argument("--schema accepts at most one run file");
            std::cout << schema(files.empty() ? std::filesystem::path{} : files.front(), base_schema);
            return 0;
        }
        if (files.empty()) throw std::invalid_argument("No TOML run file specified");
        std::vector<alps::run_configuration> runs;
        std::set<std::filesystem::path> protected_paths, destinations;
        for (auto const& file : files) protected_paths.insert(std::filesystem::weakly_canonical(file));
        for (auto const& file : files) {
            auto run = alps::load_run_configuration(file, schema(file, base_schema));
            if (run.parameters.exists("LATTICE") == run.parameters.exists("GRAPH"))
                throw std::invalid_argument("Specify exactly one parameters.LATTICE or parameters.GRAPH");
            for (auto const& [key, fallback] : {std::pair{"lattice_library", "lattices.xml"}, std::pair{"model_library", "models.xml"}})
                run.input[key] = std::filesystem::weakly_canonical(
                    alps::search_xml_library_path(run.input.value_or<std::string>(key, fallback))).string();
            for (auto const& [key, value] : run.input)
                for (auto const& path : alps::run_paths(value)) {
                    if (!std::filesystem::is_regular_file(path))
                        throw std::invalid_argument("Missing input." + key + " file: " + path);
                    protected_paths.insert(std::filesystem::weakly_canonical(path));
                }
            const auto results = std::filesystem::weakly_canonical(run.output["results"].as<std::string>());
            if (std::filesystem::is_directory(results) || !std::filesystem::is_directory(results.parent_path()))
                throw std::invalid_argument("output.results must name a file in an existing directory");
            // The lattice, model and measurement operators are checked
            // without building the Hamiltonian.
            App::with_task(run, [](auto const&) {});
            runs.push_back(std::move(run));
        }
        for (auto const& run : runs) {
            const auto path = std::filesystem::weakly_canonical(run.output["results"].as<std::string>());
            if (protected_paths.count(path) || !destinations.insert(path).second)
                throw std::invalid_argument("Output paths must be distinct and must not replace any run or input file");
        }
        if (validate) {
            for (auto const& file : files)
                std::cout << "Valid " << application << " configuration: " << file.string() << '\n';
            return 0;
        }
        App::print_copyright(std::cout);
        for (auto const& run : runs)
            App::with_task(run, [&](auto& task) {
                task.run();
                alps::hdf5::save_checkpoint(run.output["results"].as<std::string>(), [&](alps::hdf5::archive& archive) {
                    archive["/parameters"] << run.parameters;
                    archive["/run_config"] << run;
                    task.save(archive);
                });
            });
        return 0;
    } catch (std::exception const& error) {
        std::cerr << application << ": " << error.what() << '\n';
        return 1;
    }
}
} // namespace lattice_model
