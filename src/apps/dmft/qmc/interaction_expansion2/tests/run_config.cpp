// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#include <alps/ctint.hpp>
#include <alps/hdf5/vector.hpp>
#include <alps/hdf5/complex.hpp>
#include "interaction_expansion.hpp"
#include <filesystem>
#include <fstream>
#include <limits>
#ifdef ALPS_HAVE_MPI
#include <boost/mpi/environment.hpp>
#endif

void require(bool value, const char *message) {
    if (!value) throw std::runtime_error(message);
}
template <class F> void rejects(F f, const std::string &needle) {
    try { f(); }
    catch (const std::exception &error) {
        require(std::string(error.what()).find(needle) != std::string::npos, error.what());
        return;
    }
    throw std::runtime_error("invalid CT-INT configuration accepted");
}
alps::run_configuration configuration() {
    alps::run_configuration run;
    run.parameters["BETA"] = 2.0;
    run.parameters["U"] = 0.0;
    run.parameters["MU"] = 0.0;
    run.parameters["ALPHA"] = -0.01;
    run.parameters["N"] = 8;
    run.parameters["NMATSUBARA"] = 4;
    run.parameters["SWEEPS"] = 3;
    run.parameters["THERMALIZATION"] = 3;
    run.parameters["MEASUREMENT_PERIOD"] = 2;
    run.input["atomic"] = true;
    run.output["results"] = "ctint-contract-results.h5";
    run.execution["time_limit"] = 0;
    return run;
}
int main(int argc, char **argv) {
#ifdef ALPS_HAVE_MPI
    boost::mpi::environment environment(argc, argv);
#endif
    auto run = configuration();
    alps::ctint::prepare_run(run);
    require(run.application == "ctint", "schema identity missing");
    require(run.parameters["NMATSUBARA_MEASUREMENTS"].as<int>() == 4, "frequency default");
    require(run.parameters["NSELF"].as<int>() == 80, "histogram default");
    require(!run.parameters.exists("SEED") && !run.parameters.exists("MAX_TIME") &&
            !run.parameters.exists("INFILE"), "orchestration leaked into physics parameters");
    require(run.execution["seed"].as<int>() == 42, "seed default changed");
    require(run.origins.at("parameters.NSELF") == "derived", "derived provenance missing");
    for (const auto *key : {"MAX_TIME", "INFILE", "N_TAU", "N_MATSUBARA", "N_ORBITALS", "SEED", "BINNUMBER"}) {
        auto invalid = configuration();
        invalid.parameters[key] = 1;
        rejects([&] { alps::ctint::prepare_run(invalid); }, key);
    }
    auto invalid = configuration();
    invalid.parameters["BETA"] = 0.0;
    rejects([&] { alps::ctint::prepare_run(invalid); }, "BETA");
    invalid = configuration();
    invalid.parameters["NMATSUBARA_MEASUREMENTS"] = 5;
    rejects([&] { alps::ctint::prepare_run(invalid); }, "NMATSUBARA_MEASUREMENTS");
    invalid = configuration();
    invalid.parameters["MEASUREMENT_PERIOD"] = 0;
    rejects([&] { alps::ctint::prepare_run(invalid); }, "MEASUREMENT_PERIOD");
    invalid = configuration();
    invalid.parameters["FLAVORS"] = 4;
    rejects([&] { alps::ctint::prepare_run(invalid); }, "FLAVORS");
    invalid = configuration();
    invalid.input["atomic"] = false;
    rejects([&] { alps::ctint::prepare_run(invalid); }, "either input.g0");
    invalid = configuration();
    invalid.execution["check_interval"] = 0.0;
    rejects([&] { alps::ctint::prepare_run(invalid); }, "check_interval");

    // The frequency measurement limit is distinct from the input extent.
    const auto input_file = std::filesystem::absolute("ctint-contract-g0.h5").string();
    const auto write_input = [&](std::size_t count, bool nonfinite = false) {
        alps::hdf5::archive archive(input_file, "w");
        std::vector<std::complex<double>> values(count, {0.0, -0.5});
        if (nonfinite) values[0] = {std::numeric_limits<double>::infinity(), 0.0};
        archive["/G0_0"] << values;
        archive["/G0_1"] << values;
    };
    auto file_run = configuration();
    file_run.input["atomic"] = false;
    file_run.input["g0"] = input_file;
    file_run.parameters["NMATSUBARA_MEASUREMENTS"] = 2;
    write_input(4);
    alps::ctint::prepare_run(file_run);
    write_input(3);
    rejects([&] { alps::ctint::prepare_run(file_run); }, "vector of length 4");
    {
        alps::hdf5::archive archive(input_file, "w");
        const std::vector<double> malformed(12, 0.0);
        archive.write("/G0_0", malformed.data(), std::vector<std::size_t>{4, 3});
        archive.set_complex("/G0_0");
    }
    rejects([&] { alps::ctint::prepare_run(file_run); }, "complex vector");
    write_input(4, true);
    rejects([&] { alps::ctint::prepare_run(file_run); }, "must be finite");
    write_input(4);
    file_run.output["results"] = input_file;
    rejects([&] { alps::ctint::prepare_run(file_run); }, "replace input.g0");
    std::filesystem::remove(input_file);

    // Explicit thermalization must finish; prethermalization observations are excluded.
    HubbardInteractionExpansion simulation(run, 0);
    require(!simulation.is_thermalized(), "thermalized before any updates");
    simulation.update(); simulation.measure();
    require(simulation.collect_results()["Sign"].count() == 0, "warmup counted as measurement");
    simulation.update(); simulation.measure();
    require(simulation.is_thermalized() && simulation.fraction_completed() > 0.0,
            "thermalization never completed");
    simulation.update(); simulation.measure();
    require(simulation.fraction_completed() >= 1.0, "finite sweeps never complete");
    require(simulation.collect_results()["Sign"].count() == 2, "wrong measurement count");

    alps::solvers::ctint(run);
    {
        alps::hdf5::archive archive(run.output["results"].as<std::string>(), "r");
        alps::params parameters, execution;
        std::string application;
        archive["/parameters"] >> parameters;
        archive["/run_config/execution"] >> execution;
        archive["/run_config/application"] >> application;
        require(application == "ctint" && !parameters.exists("MAX_TIME"), "saved run sections");
        require(execution["seed"].as<int>() == 42, "saved execution provenance");
        std::vector<std::complex<double>> values;
        archive["/G_omega/0/mean/value"] >> values;
        require(values.size() == 4, "missing Green function output");
        require(std::abs(values[0] - std::complex<double>(0.0, -2.0 / std::acos(-1.0))) < 1e-12,
                "noninteracting Green function changed");
    }
    std::filesystem::remove(run.output["results"].as<std::string>());
    require(!std::filesystem::exists("matrix_size"), "implicit text output created");
}
