// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#include "../simulation.h"
#include <boost/math/special_functions/bessel.hpp>
#include <iostream>

namespace {
void require(bool condition, char const* message) {
    if (!condition) throw std::runtime_error(message);
}
void agree(alps::alea::batch_result<double> const& result, double expected, char const* message) {
    if (!(std::abs(result.mean()(0) - expected) <= 6. * result.stderror()(0) + 1e-3)) {
        std::cerr << message << ": measured " << result.mean()(0) << " +/- " << result.stderror()(0)
                  << ", expected " << expected << '\n';
        throw std::runtime_error(message);
    }
}
void step(simplemc::simulation& run, uint64_t count) {
    for (uint64_t i = 0; i != count; ++i) { run.update(); run.measure(); }
}
void equal(simplemc::simulation const& expected, simplemc::simulation const& actual) {
    require(expected.completed_sweeps() == actual.completed_sweeps(), "restart changed progress");
    auto reference = expected.collect_results(), resumed = actual.collect_results();
    require(reference.size() == resumed.size(), "restart changed result names");
    for (auto const& entry : reference)
        require(entry.second == resumed.at(entry.first), "restart changed native batch evidence");
}
}

int main() {
    auto directory = std::filesystem::temp_directory_path() / boost::filesystem::unique_path("simplemc-tests-%%%%-%%%%").string();
    std::filesystem::create_directory(directory);
    struct cleanup {
        std::filesystem::path path;
        ~cleanup() { std::error_code ignored; std::filesystem::remove_all(path, ignored); }
    } remove{directory};
    auto library = directory / "graphs.xml";
    {
        std::ofstream xml(library);
        xml << R"(<LATTICES>
<GRAPH name="self" vertices="1" dimension="1">
<VERTEX id="1"><COORDINATE>0</COORDINATE></VERTEX>
<EDGE source="1" target="1" type="0"/>
</GRAPH>
<GRAPH name="pair" vertices="2" dimension="1">
<VERTEX id="1"><COORDINATE>0</COORDINATE></VERTEX>
<VERTEX id="2"><COORDINATE>1</COORDINATE></VERTEX>
<EDGE source="1" target="2" type="0"/>
<EDGE source="1" target="2" type="1"/>
</GRAPH>
<GRAPH name="mixed" vertices="3" dimension="1">
<VERTEX id="1"><COORDINATE>0</COORDINATE></VERTEX>
<VERTEX id="2"><COORDINATE>1</COORDINATE></VERTEX>
<VERTEX id="3"><COORDINATE>2</COORDINATE></VERTEX>
<EDGE source="1" target="2" type="0"/>
<EDGE source="2" target="3" type="1"/>
<EDGE source="1" target="3" type="2"/>
<EDGE source="3" target="3" type="3"/>
</GRAPH>
</LATTICES>)";
    }
    alps::params p;
    p["LATTICE_LIBRARY"] = library.string(); p["GRAPH"] = "self";
    p["THERMALIZATION"] = 500; p["SWEEPS"] = 60000; p["SEED"] = 42;
    p["T"] = .8; p["H"] = .9; p["J"] = 1.7;
    for (auto const& model : {"ising", "xy", "heisenberg"}) {
        p["ALGORITHM"] = model;
        auto changed = p; changed["J"] = 3.4;
        simplemc::simulation reference(p, 16), changed_loop(changed, 16);
        step(reference, 757); step(changed_loop, 757);
        auto projection = std::string("Magnetization Density") +
            (std::string(model) == "ising" ? "" : std::string(model) == "xy" ? " X" : " Z");
        require(reference.collect_results().at(projection) == changed_loop.collect_results().at(projection)
                && reference.get_random()() == changed_loop.get_random()(), "self-loop J changed the spin trajectory");
        simplemc::simulation run(p, 64);
        step(run, 60500);
        auto results = run.collect_results();
        double x = .9 / .8;
        double magnetization = std::string(model) == "ising" ? std::tanh(x)
            : std::string(model) == "xy" ? boost::math::cyl_bessel_i(1, x) / boost::math::cyl_bessel_i(0, x)
            : 1. / std::tanh(x) - 1. / x;
        agree(results.at(projection), magnetization, "self-loop changed the field Boltzmann distribution");
        agree(results.at("Energy"), -1.7 - .9 * magnetization, "self-loop energy or field axis is incorrect");
        require(run.measurement_count() == 60000 && run.fraction_completed() == 1., "incorrect thermalization count");
        auto snapshot = directory / (std::string(model) + ".vtk");
        run.snapshot(snapshot);
        auto size = std::filesystem::file_size(snapshot);
        bool rejected = false;
        try { run.snapshot(snapshot); } catch (std::filesystem::filesystem_error const&) { rejected = true; }
        require(rejected && std::filesystem::file_size(snapshot) == size, "snapshot overwrote an existing file");
    }

    // Independent analytic pair oracles check the dot-product energy and
    // magnetization normalization for both continuous-spin models.
    p["GRAPH"] = "pair"; p["H"] = 0.; p["T"] = .7; p["J0"] = .6; p["J1"] = -.1;
    for (auto const& model : {"xy", "heisenberg"}) {
        p["ALGORITHM"] = model;
        simplemc::simulation run(p, 64);
        step(run, 60500);
        double x = .5 / .7;
        double correlation = std::string(model) == "xy"
            ? boost::math::cyl_bessel_i(1, x) / boost::math::cyl_bessel_i(0, x)
            : 1. / std::tanh(x) - 1. / x;
        auto results = run.collect_results();
        agree(results.at("Energy"), -.5 * correlation, "continuous-spin pair energy disagrees with exact integral");
        agree(results.at("Magnetization Density^2"), (1 + correlation) / 2., "continuous-spin pair magnetization is incorrect");
    }

    // Enumerating all Ising states is independent of the engine's update rule.
    p["ALGORITHM"] = "ising"; p["GRAPH"] = "mixed"; p["T"] = 1.2; p["H"] = .45;
    p["J0"] = .7; p["J1"] = -.4; p["J2"] = .3; p["J3"] = 1.1;
    double partition = 0, e = 0, e2 = 0, m2 = 0, m4 = 0;
    for (int state = 0; state != 8; ++state) {
        double a = state & 1 ? 1 : -1, b = state & 2 ? 1 : -1, c = state & 4 ? 1 : -1;
        double energy = -.7*a*b + .4*b*c - .3*a*c - 1.1 - .45*(a+b+c);
        double m = (a+b+c)/3., weight = std::exp(-energy/1.2);
        partition += weight; e += weight*energy; e2 += weight*energy*energy;
        m2 += weight*m*m; m4 += weight*m*m*m*m;
    }
    simplemc::simulation ising(p, 64);
    step(ising, 60500);
    auto results = simplemc::simulation::derive(ising.collect_results(), p);
    agree(results.at("Energy"), e/partition, "Ising mixed-bond energy disagrees with enumeration");
    agree(results.at("Energy^2"), e2/partition, "Ising energy second moment disagrees with enumeration");
    agree(results.at("Specific Heat"), (e2/partition - std::pow(e/partition, 2))/(3*1.2*1.2), "specific heat lost joint moment covariance");
    agree(results.at("Binder Ratio of Magnetization"), m2*m2/(partition*m4), "Binder ratio disagrees with enumeration");

    p["THERMALIZATION"] = 7; p["SWEEPS"] = 1037;
    for (auto const& model : {"ising", "xy", "heisenberg"}) {
        p["ALGORITHM"] = model;
        simplemc::simulation uninterrupted(p, 16, 3), stopped(p, 16, 3), resumed(p, 16, 3);
        step(uninterrupted, 1044); step(stopped, 148);
        auto checkpoint = (directory / "state.h5").string();
        stopped.save(checkpoint);
        resumed.load(checkpoint);
        step(resumed, 1044 - 148);
        equal(uninterrupted, resumed);
        require(uninterrupted.get_random()() == resumed.get_random()(), "restart changed RNG continuation");
        {
            alps::hdf5::archive ar(checkpoint, "a");
            ar["/simulation/realizations/0/clones/0/checkpoint/chain_id"] << uint64_t(4);
        }
        bool rejected = false;
        try { resumed.load(checkpoint); } catch (std::invalid_argument const&) { rejected = true; }
        require(rejected, "checkpoint accepted a different chain ID");
        equal(uninterrupted, resumed);
        stopped.save(checkpoint);
        {
            alps::hdf5::archive ar(checkpoint, "a");
            ar["/simulation/realizations/0/clones/0/checkpoint/spins"] << std::vector<std::array<double,3>>(3, {0.,0.,0.});
        }
        rejected = false;
        try { resumed.load(checkpoint); } catch (std::invalid_argument const&) { rejected = true; }
        require(rejected, "checkpoint accepted a non-unit spin");
        equal(uninterrupted, resumed);
        stopped.save(checkpoint);
        simplemc::simulation wrong_bins(p, 32, 3);
        rejected = false;
        try { wrong_bins.load(checkpoint); } catch (std::invalid_argument const&) { rejected = true; }
        require(rejected && wrong_bins.completed_sweeps() == 0, "checkpoint accepted a different batching configuration");
        stopped.save(checkpoint);
        {
            alps::hdf5::archive ar(checkpoint, "a");
            auto path = std::string("/simulation/realizations/0/clones/0/measurements");
            for (auto const& name : ar.list_children(path)) {
                uint64_t position;
                ar[path + "/" + name + "/cursor/level_position"] >> position;
                ar[path + "/" + name + "/cursor/level_position"] << (position + 1) % 8;
            }
        }
        rejected = false;
        try { resumed.load(checkpoint); } catch (std::runtime_error const&) { rejected = true; }
        require(rejected, "checkpoint accepted aligned but impossible moment cursors");
        equal(uninterrupted, resumed);
    }

    // No warm-up still records the first sweep. Unequal chains retain their
    // own bins and sample weights, including a partially filled last bin.
    p["THERMALIZATION"] = 0; p["SWEEPS"] = 29;
    simplemc::simulation first(p, 8, 0), second(p, 8, 1);
    step(first, 1); require(first.measurement_count() == 1, "zero warm-up lost the first sample");
    step(first, 28); step(second, 17);
    auto a = first.collect_results(), b = second.collect_results();
    auto merged = simplemc::simulation::merge({a, b});
    for (auto const& entry : merged) {
        auto const& data = entry.second;
        require(data.count() == 46 && data.num_batches() == 16, "merge discarded independent chain evidence");
        require(data.store().batch().leftCols(8) == a.at(entry.first).store().batch()
                && data.store().batch().rightCols(8) == b.at(entry.first).store().batch(), "merge added independent chain bins");
        require(std::abs(data.mean()(0) - (29*a.at(entry.first).mean()(0) + 17*b.at(entry.first).mean()(0))/46.) < 1e-12,
                "merge weighted unequal chains incorrectly");
    }
    auto excessive = p; excessive["J0"] = 1e150;
    bool rejected = false;
    try { simplemc::simulation invalid(excessive); } catch (std::invalid_argument const&) { rejected = true; }
    require(rejected, "preparation accepted overflowing energy-squared uncertainty");
    auto frozen = p; frozen["ALGORITHM"] = "ising"; frozen["GRAPH"] = "self";
    frozen["T"] = 1e-200; frozen["H"] = 1.; frozen["J0"] = 0.;
    simplemc::simulation cold(frozen, 8); step(cold, 29);
    auto cold_results = simplemc::simulation::derive(cold.collect_results(), frozen);
    require(cold_results.at("Specific Heat").mean().allFinite(), "zero heat capacity overflowed at finite inverse temperature");
    std::cout << "simplemc native physics, restart, covariance and chain contracts passed\n";
}
