/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2013 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


#include "ising.hpp"

#include <alps/mc/driver.hpp>

#include <iostream>
#include <stdexcept>
#include <vector>

// alps::mc::main runs execution.chains independent chains of each TOML run
// file, spread over the MPI ranks, and checkpoints and restarts them:
//
//     mpiexec -n 2 ./mpi param.0.toml param.1.toml
//
// Several run files form a parameter scan; they run one after another, each
// on all ranks. `./mpi --schema` prints the schema below.
constexpr char schema[] = R"toml(
application = "ising-tutorial"
schema_version = 1
[parameters.L]
type = "int64"
required = true
min = 1
[parameters.T]
type = "float64"
required = true
min = 0.0
[parameters.THERMALIZATION]
type = "int64"
required = true
min = 0
[parameters.SWEEPS]
type = "int64"
required = true
min = 1
[input.checkpoint]
type = "path"
[output.results]
type = "path"
required = true
[output.checkpoint]
type = "path"
[execution.chains]
type = "int64"
default = 1
min = 1
max = 2147483647
[execution.seed]
type = "int64"
default = 42
min = 0
max = 2147483647
[execution.rng]
type = "string"
default = "mt19937"
choices = ["mt19937", "lagged_fibonacci607"]
[execution.bins]
type = "int64"
default = 64
min = 2
[execution.time_limit]
type = "float64"
default = 0.0
min = 0.0
[execution.checkpoint_interval]
type = "float64"
default = 3600.0
min = 0.0
[execution.max_sweeps]
type = "int64"
default = 0
min = 0
)toml";

int main(int argc, char *argv[]) {
    auto prepare = [](alps::params & parameters, alps::run_configuration const &) {
        if (!(parameters["T"].as<double>() > 0))
            throw std::invalid_argument("T must be positive");
    };
    // Rank 0 publishes once all chains of a run are complete or stopped.
    auto publish = [](alps::run_configuration const & run, auto const & chains, alps::params const &) {
        std::vector<alps::mc::batch_results> results;
        for (auto const & chain : chains)
            results.push_back(chain->collect_results());
        auto const pooled = alps::mc::pool(results);
        for (auto const & entry : pooled)
            std::cout << entry.first << ": " << entry.second << '\n';
        alps::hdf5::save_checkpoint(run.output["results"].as<std::string>(), [&](alps::hdf5::archive & archive) {
            alps::save_results(pooled, run.parameters, archive, "/simulation/results");
            archive["/run_config"] << run;
        });
    };
    return alps::mc::main<ising_sim>(argc, argv, "mpi", schema, {}, prepare, publish);
}
