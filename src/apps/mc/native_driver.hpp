// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#pragma once
#include "measurements.hpp"
#include "parallel.hpp"
#include <alps/alea/convert.hpp>
#include <alps/alea/batch.hpp>
#include <alps/ngs/signal.hpp>
#include <alps/ngs/api.hpp>
#include <alps/run_config.hpp>
#include <alps/parser/xslt_path.h>
#include <toml++/toml.hpp>
#include <chrono>
#include <climits>
#include <filesystem>
#include <functional>
#include <iostream>
#include <limits>
#include <memory>
#include <set>
#include <sstream>

// Private application orchestration: retain independent, weighted batches
// before any nonlinear estimator joins the aligned physical moments.
namespace native_mc {
using batch_results = std::map<std::string,alps::alea::batch_result<double>>;
inline batch_results pool(std::vector<batch_results> const& chains) {
    if (chains.empty()) return {};
    for (auto const& chain : chains) {
        if (chain.empty() || chain.size() != chains.front().size())
            throw std::invalid_argument("Inconsistent chain result names");
        auto const& reference = chain.begin()->second;
        for (auto const& [name, first] : chains.front()) {
            auto const& result = chain.at(name);
            if (!result.valid() || !result.size() || result.size() != first.size()
                    || result.num_batches() != reference.num_batches()
                    || result.store().count() != reference.store().count() || !result.store().batch().allFinite())
                throw std::invalid_argument("Chain moments are not aligned");
            for (Eigen::Index i = 0; i < result.store().count().size(); ++i)
                if (!result.store().count()(i) && !result.store().batch().col(i).isZero(0))
                    throw std::invalid_argument("Empty chain batch contains observations");
        }
    }
    batch_results merged;
    for (auto const& [name, first] : chains.front()) {
        std::vector<alps::alea::batch_result<double>> runs;
        for (auto const& chain:chains) runs.push_back(chain.at(name));
        merged.emplace(name,alps::alea::merge(runs));
    }
    return merged;
}

using parameter_type = std::function<char const*(std::string const&, toml::node const&)>;
inline std::string schema(std::filesystem::path const& file, char const* base, parameter_type const& type) {
    auto definitions = toml::parse(base);
    if (!file.empty()) {
        const auto document = toml::parse_file(file.string());
        if (const auto* parameters = document["parameters"].as_table())
            for (auto const& [name, node] : *parameters) {
                const auto key = std::string(name.str());
                if (definitions["parameters"].as_table()->contains(key)) continue;
                if (key == "NUM_CLONES" || key == "SEED" || key == "RNG" || key == "SNAPSHOT_INTERVAL" ||
                    key == "LATTICE_LIBRARY" || key == "WORKER" || key == "WORKER_SEED" ||
                    key == "DISORDER_SEED" || key == "DISORDERSEED" || key == "ERROR_VARIABLE" ||
                    key == "ERROR_LIMIT" || key == "PRINT_SWEEPS")
                    throw std::invalid_argument("Retired parameter " + key + "; use typed input/output/execution fields");
                const char* value_type = type(key, node);
                if (!value_type) value_type = node.is_integer() ? "int64" : node.is_floating_point() ?
                    "float64" : node.is_string() ? "string" : node.is_boolean() ? "bool" : nullptr;
                if (!value_type) throw std::invalid_argument("Graph parameter " + key + " must be a scalar");
                definitions["parameters"].as_table()->insert(key, toml::table{{"type", value_type}});
            }
    }
    std::ostringstream output;
    output << definitions;
    return output.str();
}

template<class Prepare>
alps::params parameters(alps::run_configuration const& run, Prepare const& prepare) {
    auto p = run.parameters;
    p["SEED"] = run.execution["seed"];
    p["RNG"] = run.execution["rng"];
    p["DISORDER_SEED"] = run.execution.value_or<std::uint64_t>("disorder_seed", run.execution["seed"].as<std::uint64_t>());
    p["LATTICE_LIBRARY"] = run.input["lattice_library"];
    prepare(p, run);
    return p;
}

template<class Simulation> using chains_type = std::vector<std::unique_ptr<Simulation>>;
template<class Simulation, class Prepare>
chains_type<Simulation> prepare_chains(alps::run_configuration const& run, Prepare const& prepare, parallel const& group) {
    auto const& execution = run.execution;
    const auto chains = execution["chains"].as<std::size_t>();
    const auto seed = execution["seed"].as<std::uint64_t>();
    if (chains - 1 > INT_MAX - seed)
        throw std::invalid_argument("execution.seed + chains - 1 exceeds the supported RNG seed range");
    if (execution["bins"].as<std::size_t>() % 2)
        throw std::invalid_argument("execution.bins must be even and at least two");
    if (run.parameters.exists("LATTICE") == run.parameters.exists("GRAPH"))
        throw std::invalid_argument("Specify exactly one parameters.LATTICE or parameters.GRAPH");
    if (execution.value_or<std::uint64_t>("snapshot_interval", 0) && !run.output.exists("snapshot_prefix"))
        throw std::invalid_argument("Snapshots require output.snapshot_prefix");
    for (auto const& [key, value] : run.input)
        if (!std::filesystem::is_regular_file(value.as<std::string>()))
            throw std::invalid_argument("Missing input." + key + " file: " + value.as<std::string>());
    for (auto const& [key, value] : run.output) {
        const auto path = std::filesystem::weakly_canonical(value.as<std::string>());
        if (std::filesystem::is_directory(path) || !std::filesystem::is_directory(path.parent_path()))
            throw std::invalid_argument("output." + key + " must name a file in an existing directory");
    }
    const alps::params p = parameters(run, prepare);
    if (p["SWEEPS"].as<std::uint64_t>() > std::numeric_limits<std::uint64_t>::max() - p["THERMALIZATION"].as<std::uint64_t>())
        throw std::invalid_argument("Total sweep count exceeds the counter range");
    chains_type<Simulation> simulations;
    for (std::size_t id = 0; id < chains; ++id)
        simulations.push_back(group.rank()==0 || group.owns(id)
            ? std::make_unique<Simulation>(p, execution["bins"].as<std::size_t>(), id) : nullptr);
    if (run.input.exists("checkpoint")) {
        alps::hdf5::archive archive(run.input["checkpoint"].as<std::string>());
        alps::run_configuration saved;
        archive["/run_config"] >> saved;
        if (saved.application != run.application || saved.schema_version != run.schema_version ||
            Simulation::checkpoint_parameters(parameters(saved, prepare)) != Simulation::checkpoint_parameters(p)
            || saved.execution["chains"].as<std::size_t>() != chains ||
            saved.execution["bins"].as<std::size_t>() != execution["bins"].as<std::size_t>())
            throw std::invalid_argument("Checkpoint model, graph, seed, chains and bins must match this run");
        const auto clones = "/simulation/realizations/0/clones/";
        if (archive.list_children(clones).size() != chains)
            throw std::invalid_argument("Checkpoint clone set does not match execution.chains");
        for (std::size_t id = 0; id < chains; ++id)
            if (!archive.is_group(clones + std::to_string(id)))
                throw std::invalid_argument("Unexpected checkpoint clone name or type");
        for (std::size_t id = 0; id < chains; ++id) if (simulations[id]) {
            archive.set_context(clones + std::to_string(id));
            simulations[id]->load(archive);
        }
    }
    return simulations;
}

template<class Simulation>
void checkpoint(alps::run_configuration const& run, chains_type<Simulation> const& chains) {
    if (!run.output.exists("checkpoint")) return;
    alps::hdf5::save_checkpoint(run.output["checkpoint"].as<std::string>(), [&](alps::hdf5::archive& archive) {
        archive["/run_config"] << run;
        for (std::size_t id = 0; id < chains.size(); ++id) {
            archive.set_context("/simulation/realizations/0/clones/" + std::to_string(id));
            chains[id]->save(archive);
        }
    });
}

template<class Simulation> using snapshot_type = std::function<void(Simulation const&, std::filesystem::path const&)>;
template<class Simulation, class Prepare, class Derive>
void execute(alps::run_configuration const& run, chains_type<Simulation>& chains,
             Prepare const& prepare, Derive const& derive, snapshot_type<Simulation> const& snapshot, parallel const& group) {
    using clock = std::chrono::steady_clock;
    const auto started = clock::now();
    auto last_checkpoint = started;
    const auto limit = run.execution["time_limit"].as<double>();
    const auto interval = run.execution["checkpoint_interval"].as<double>();
    const auto budget = run.execution["max_sweeps"].as<std::uint64_t>();
    const auto snapshots = run.execution.value_or<std::uint64_t>("snapshot_interval", 0);
    std::vector<std::uint64_t> initial;
    for (auto const& chain : chains) initial.push_back(chain ? chain->completed_sweeps() : 0);
    alps::ngs::signal signal;
    auto stopped=[&] { return !signal.empty() || (limit && std::chrono::duration<double>(clock::now()-started).count()>=limit); };
    bool active=true;
    while (active && !group.any(stopped())) {
        group.checked([&] {
            // Amortize communication without delaying local stopping checks.
            for (int step=0;step<32 && !stopped();++step) {
                active=false;
                for (std::size_t id=0;id<chains.size();++id) if (group.owns(id)) {
                    auto& chain=*chains[id];
                    if (chain.fraction_completed()>=1. || (budget && chain.completed_sweeps()-initial[id]>=budget)) continue;
                    active=true;
                    chain.update(); chain.measure();
                    if (snapshots && chain.completed_sweeps()%snapshots==0)
                        snapshot(chain,run.output["snapshot_prefix"].as<std::string>()+".clone"+
                            std::to_string(id+1)+"."+std::to_string(chain.completed_sweeps())+".vtk");
                }
                if (!active) break;
            }
        });
        active=group.any(active);
        if (group.any(interval && run.output.exists("checkpoint") &&
                std::chrono::duration<double>(clock::now()-last_checkpoint).count()>=interval)) {
            group.synchronize(chains);
            group.checked([&] { if (group.rank()==0) checkpoint(run,chains); });
            last_checkpoint=clock::now();
        }
    }
    group.synchronize(chains);
    group.checked([&] {
    if (group.rank()!=0) return;
    checkpoint(run,chains);
    std::vector<batch_results> raw;
    for (auto const& chain : chains) raw.push_back(chain->collect_results());
    const auto results = derive(pool(raw), parameters(run, prepare));
    alps::hdf5::save_checkpoint(run.output["results"].as<std::string>(), [&](alps::hdf5::archive& archive) {
        alps::save_results(results, run.parameters, archive, "/simulation/results");
        archive["/run_config"] << run;
        for (std::size_t id = 0; id < chains.size(); ++id) {
            archive["/simulation/realizations/0/clones/" + std::to_string(id) + "/completed_sweeps"] << chains[id]->completed_sweeps();
            archive["/simulation/realizations/0/clones/" + std::to_string(id) + "/measurements"] << chains[id]->measurement_count();
            save_diagnostics(*chains[id],archive,"/simulation/realizations/0/clones/"+std::to_string(id),
                             native_mc::batch_names(chains[id]->get_measurements()));
        }
    });
    });
}

template<class Simulation, class Prepare, class Derive>
int main(int argc, char** argv, char const* application, char const* base_schema,
         parameter_type type, Prepare prepare, Derive derive, snapshot_type<Simulation> snapshot = {}) {
    try {
        bool validate = false, show_schema = false;
        std::vector<std::filesystem::path> files;
        for (int i = 1; i < argc; ++i) {
            const std::string argument(argv[i]);
            if (argument == "--help" || argument == "-h") {
                std::cout << "Usage: " << application << " [--validate] run.toml [run.toml ...] | " << application << " --schema [run.toml]\n"
                    "Independent chains, RNG seed and stopping belong in [execution].\n"
                    "Use input.lattice_library/checkpoint and output.results/checkpoint.\n";
                return 0;
            }
            if (argument == "--validate") validate = true;
            else if (argument == "--schema") show_schema = true;
            else if (argument.empty() || argument.front() == '-') throw std::invalid_argument("Unknown option: " + argument);
            else files.emplace_back(argument);
        }
        if (show_schema) {
            if (files.size() > 1 || validate) throw std::invalid_argument("--schema accepts at most one run file");
            std::cout << schema(files.empty() ? std::filesystem::path{} : files.front(), base_schema, type);
            return 0;
        }
        if (files.empty()) throw std::invalid_argument("No TOML run file specified");
#ifdef ALPS_HAVE_MPI
        boost::mpi::environment environment(argc, argv, false);
#endif
        parallel group;
        std::vector<alps::run_configuration> runs;
        std::vector<chains_type<Simulation>> simulations;
        std::set<std::filesystem::path> protected_paths, destinations;
        group.checked([&] {
        for (auto const& file : files) protected_paths.insert(std::filesystem::weakly_canonical(file));
        for (auto const& file : files) {
            runs.push_back(alps::load_run_configuration(file, schema(file, base_schema, type)));
            runs.back().input["lattice_library"] = std::filesystem::weakly_canonical(
                alps::search_xml_library_path(runs.back().input.value_or<std::string>("lattice_library", "lattices.xml"))).string();
            simulations.push_back(prepare_chains<Simulation>(runs.back(), prepare, group));
            for (auto const& [key, value] : runs.back().input) protected_paths.insert(std::filesystem::weakly_canonical(value.as<std::string>()));
        }
        for (auto const& run : runs)
            for (auto const& [key, value] : run.output) {
                const auto path = std::filesystem::weakly_canonical(value.as<std::string>());
                if (protected_paths.count(path) || !destinations.insert(path).second)
                    throw std::invalid_argument("Output paths must be distinct and must not replace any run or input file");
            }
        for (auto const& run : runs)
            if (run.output.exists("snapshot_prefix")) {
                const auto prefix = std::filesystem::weakly_canonical(run.output["snapshot_prefix"].as<std::string>()).string() + ".clone";
                for (auto const& path : protected_paths)
                    if (path.string().compare(0, prefix.size(), prefix) == 0)
                        throw std::invalid_argument("Snapshot prefix would replace a run or input file");
                for (auto const& path : destinations)
                    if (path.string().compare(0, prefix.size(), prefix) == 0)
                        throw std::invalid_argument("Snapshot prefix overlaps another output path or prefix");
            }
        });
        for (std::size_t index = 0; index < runs.size(); ++index) {
            if (validate) { if (group.rank()==0) std::cout << "Valid " << application << " configuration: " << files[index].string() << '\n'; }
            else execute(runs[index], simulations[index], prepare, derive, snapshot, group);
        }
        return 0;
    } catch (std::exception const& error) {
        std::cerr << application << ": " << error.what() << '\n';
        return 1;
    }
}
} // namespace native_mc
