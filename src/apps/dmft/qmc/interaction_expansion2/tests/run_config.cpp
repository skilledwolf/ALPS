// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#include <alps/ctint.hpp>
#include <alps/hdf5/vector.hpp>
#include <alps/hdf5/complex.hpp>
#include <alps/alea/hdf5.hpp>
#include "interaction_expansion.hpp"
#include "run_config.hpp"
#include <filesystem>
#include <fstream>
#include <limits>
#include <array>
#include <cmath>
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
// Exercise the actual production registry without changing the Markov chain.
class measured_simulation : public InteractionExpansion {
public:
    using InteractionExpansion::InteractionExpansion;
    void density_sample(double sample_sign, std::array<double, 2> const& density) {
        sign = sample_sign;
        record_measurement("densities", std::valarray<double>{sign * density[0], sign * density[1]});
        record_measurement("Sign", sign);
        record_measurement("VertexInsertion", 1.);
    }
    auto raw_density() const { return measurements.at("densities").accumulator.result(); }
    void rebuild_negative_configuration() {
        // U=-4, beta=2 and atomic G0=-1/2 give a first-vertex weight -2.0808,
        // independently of the sampled vertex time and symmetrized alpha.
        auto proposal = try_add();
        require(std::abs(proposal + 2.0808) < 1e-12, "negative insertion fixture changed");
        require(proposal < -1., "negative vertex is not unconditionally accepted");
        perform_add();
        sign *= proposal < 0. ? -1. : 1.;
        auto before = M;
        require(vertices.size() == 1 && sign == -1., "negative configuration was not accepted");
        reset_perturbation_series();
        require(vertices.size() == 1 && sign == -1., "matrix rebuild changed the accepted sign");
        for (std::size_t flavor = 0; flavor < M.size(); ++flavor) {
            require(num_rows(M[flavor].matrix()) == 1, "matrix rebuild lost the vertex");
            require(std::abs(M[flavor].matrix()(0, 0) - before[flavor].matrix()(0, 0)) < 1e-12,
                    "matrix rebuild changed the inverse propagator");
        }
    }
};

void statistics_contract(alps::run_configuration run) {
    run.execution["bins"] = 8;
    measured_simulation simulation(run, 0);
    auto empty = simulation.collect_results();
    require(empty.at("densities").count() == 0 && empty.at("densities").size() == 2,
            "empty signed result lost its physical shape");
    require(empty.at("Wk_real_0_0_0").size() == 4, "frequency shape includes sign");
    simulation.density_sample(-1., {2.5, -7.});
    auto single = simulation.collect_results().at("densities");
    require(single.count() == 1 && single.mean()(0) == 2.5 && single.mean()(1) == -7.,
            "single-bin signed mean was discarded");
    require(!single.stderror().array().isFinite().any(), "single-bin error claims independent information");
    require(simulation.collect_results().at("VertexInsertion").mean()(0) == 1.,
            "unsigned diagnostics were sign-normalized");

    measured_simulation sampled(run, 0);
    std::array<long double, 3> sum{};
    constexpr std::uint64_t samples = 173;
    for (std::uint64_t i = 0; i < samples; ++i) {
        auto sign = i % 5 == 0 ? -1. : 1.;
        std::array<double, 2> density{2. + (i % 9) * .125, -1. + (i % 7) * .5};
        sampled.density_sample(sign, density);
        for (std::size_t j = 0; j < 2; ++j) sum[j] += sign * density[j];
        sum[2] += sign;
    }
    auto raw = sampled.raw_density();
    require(raw.count() == samples, "partial bins lost samples");
    for (std::size_t j = 0; j < 3; ++j)
        require(std::abs(raw.store().batch().row(j).sum() - sum[j]) < 1e-12,
                "raw joint batch differs from independent samples");
    auto result = sampled.collect_results().at("densities");
    require(result.count() == samples && result.store().count() == raw.store().count(),
            "ratio analysis lost bin weights");
    require(sampled.raw_density() == raw, "result collection consumed live measurements");
    require(raw.store().count().minCoeff() != raw.store().count().maxCoeff(),
            "fixture did not cover unequal partial bins");
    // Independent weighted jackknife oracle, evaluated in extended precision.
    long double count2 = 0;
    for (auto weight : raw.store().count()) count2 += static_cast<long double>(weight) * weight;
    for (std::size_t j = 0; j < 2; ++j) {
        std::vector<long double> pseudo(raw.num_batches());
        long double mean = 0;
        for (std::size_t i = 0; i < raw.num_batches(); ++i) {
            auto weight = raw.store().count()(i);
            if (!weight) continue;
            auto leave_ratio = (sum[j] - raw.store().batch()(j, i))
                             / (sum[2] - raw.store().batch()(2, i));
            pseudo[i] = (samples * (sum[j] / sum[2]) - (samples - weight) * leave_ratio) / weight;
            mean += weight * pseudo[i] / samples;
        }
        long double squared = 0;
        for (std::size_t i = 0; i < pseudo.size(); ++i) {
            auto difference = pseudo[i] - mean;
            squared += raw.store().count()(i) * difference * difference;
        }
        auto variance = squared / (samples - count2 / samples);
        auto error = std::sqrt(variance * count2 / (samples * samples));
        require(std::abs(result.mean()(j) - mean) < 1e-12, "signed mean differs from weighted-bin oracle");
        require(std::abs(result.stderror()(j) - error) < 1e-12, "signed error loses joint covariance");
    }
    measured_simulation cancelled(run, 0);
    cancelled.density_sample(1., {2., 3.});
    cancelled.density_sample(-1., {2., 3.});
    auto before = cancelled.raw_density();
    rejects([&] { cancelled.collect_results(); }, "zero average sign");
    require(cancelled.raw_density() == before, "failed analysis mutated live state");
    run.parameters["HISTOGRAM_MEASUREMENT"] = true;
    measured_simulation histogram(run, 0);
    auto histogram_empty = histogram.collect_results();
    require(histogram_empty.at("W_0_0_0").size() == 81 && histogram_empty.at("Sz_0").size() == 1,
            "histogram result shapes changed");
}

void atomic_moment_contract(alps::run_configuration run) {
    const auto origins = run.origins;
    alps::ctint::prepare_run(run);
    require(run.origins == origins, "repeated preparation changed atomic moment provenance");
    for (const auto *key : {"EPS_0", "EPS_1", "EPSSQ_0", "EPSSQ_1"}) {
        require(run.parameters[key].as<double>() == 0. &&
                run.origins.at(std::string("parameters.") + key) == "derived",
                "atomic moments do not describe the zero-energy level");
        auto explicit_run = configuration();
        explicit_run.parameters[key] = 1.;
        rejects([&] { alps::ctint::prepare_run(explicit_run); }, key);
        explicit_run.parameters[key] = 0.;
        alps::ctint::prepare_run(explicit_run);
        require(explicit_run.origins.at(std::string("parameters.") + key) == "input",
                "explicit zero moment provenance was discarded");
    }
    const char *filename = "ctint-contract-run.h5";
    {
        alps::hdf5::archive archive(filename, "w");
        archive["/run"] << run;
    }
    alps::run_configuration restored;
    {
        alps::hdf5::archive archive(filename, "r");
        archive["/run"] >> restored;
    }
    alps::ctint::prepare_run(restored);
    require(restored.parameters == run.parameters && restored.origins == run.origins,
            "archived atomic moments or provenance changed on revalidation");
    std::filesystem::remove(filename);
}

void multiband_input_contract() {
    auto run = configuration();
    run.parameters["FLAVORS"] = 4;
    run.parameters["J"] = .1;
    run.parameters["U'"] = .2;
    const auto schema = alps::ctint::schema(run.parameters);
    auto resolved = alps::resolve_run_configuration(run, schema);
    require(resolved.parameters["EPS_3"].as<double>() == 0. &&
            resolved.parameters["EPSSQ_3"].as<double>() == 1., "multiband moment defaults missing");
    alps::ctint::prepare_run(run);
    for (unsigned flavor = 0; flavor < 4; ++flavor)
        for (const auto *prefix : {"EPS_", "EPSSQ_"}) {
            const auto key = std::string(prefix) + std::to_string(flavor);
            require(run.parameters[key].as<double>() == 0. && run.origins.at("parameters." + key) == "derived",
                    "multiband atomic moments were not derived");
        }
    auto invalid = configuration();
    invalid.parameters.erase("U");
    rejects([&] { alps::ctint::prepare_run(invalid); }, "parameters.U or input.interaction_matrix");
    invalid = configuration();
    invalid.parameters["FLAVORS"] = 3;
    rejects([&] { alps::ctint::prepare_run(invalid); }, "paired FLAVORS");
    invalid = configuration();
    invalid.parameters["FLAVORS"] = 4;
    invalid.parameters["EPS_4"] = 0.;
    rejects([&] { alps::ctint::prepare_run(invalid); }, "EPS_4");
    invalid = configuration();
    invalid.parameters["FLAVORS"] = 128;
    invalid.parameters["N"] = 214748363;
    rejects([&] { alps::ctint::prepare_run(invalid); }, "storage range");
    invalid = configuration();
    invalid.parameters["BETA"] = std::numeric_limits<double>::denorm_min();
    rejects([&] { alps::ctint::prepare_run(invalid); }, "Fourier arithmetic range");
    invalid = configuration();
    invalid.parameters["MU"] = std::numeric_limits<double>::max();
    invalid.parameters["H"] = std::numeric_limits<double>::max();
    rejects([&] { alps::ctint::prepare_run(invalid); }, "Fourier arithmetic range");

    const auto matrix = std::filesystem::absolute("ctint-contract-interaction.dat");
    const auto write_matrix = [&](const std::string &text) { std::ofstream(matrix) << text; };
    run = configuration();
    run.parameters["FLAVORS"] = 3;
    run.parameters.erase("U");
    run.input["interaction_matrix"] = matrix.string();
    write_matrix("0 1 1.5\n1 0 1.5\n");
    alps::ctint::prepare_run(run);  // An isolated flavor is valid.
    require(!run.parameters.exists("U"), "explicit matrix inserted an unused scalar U");
    write_matrix("");
    alps::ctint::prepare_run(run);  // A zero interaction matrix is valid.
    write_matrix("0 1 1.5\n");
    rejects([&] { alps::ctint::prepare_run(run); }, "symmetric");
    write_matrix("0 0 1.5\n");
    rejects([&] { alps::ctint::prepare_run(run); }, "zero diagonal");
    write_matrix("3 0 1.5\n");
    rejects([&] { alps::ctint::prepare_run(run); }, "Invalid index or value");
    write_matrix("0 1 nope\n");
    rejects([&] { alps::ctint::prepare_run(run); }, "Malformed input.interaction_matrix");
    write_matrix("0 1 1e200\n1 0 1e200\n");
    rejects([&] { alps::ctint::prepare_run(run); }, "interaction moments");
    write_matrix("0 1 1.5\n1 0 1.5\n");
    run.output["results"] = matrix.string();
    rejects([&] { alps::ctint::prepare_run(run); }, "replace input.interaction_matrix");
    std::filesystem::remove(matrix);
}

int main(int argc, char **argv) {
#ifdef ALPS_HAVE_MPI
    boost::mpi::environment environment(argc, argv);
#endif
    auto run = configuration();
    alps::ctint::prepare_run(run);
    atomic_moment_contract(run);
    multiband_input_contract();
    statistics_contract(run);
    auto negative_run = run;
    negative_run.parameters["U"] = -4.0;
    measured_simulation negative_configuration(negative_run, 0);
    negative_configuration.rebuild_negative_configuration();
    for (auto bins : {0, 1, 3, 127}) {
        auto invalid_bins = configuration();
        invalid_bins.execution["bins"] = bins;
        rejects([&] { alps::ctint::prepare_run(invalid_bins); }, "bins");
    }
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
    invalid.parameters["FLAVORS"] = 1;
    rejects([&] { alps::ctint::prepare_run(invalid); }, "FLAVORS");
    invalid = configuration();
    invalid.input["atomic"] = false;
    rejects([&] { alps::ctint::prepare_run(invalid); }, "either input.g0");
    invalid = configuration();
    invalid.execution["check_interval"] = 0.0;
    rejects([&] { alps::ctint::prepare_run(invalid); }, "check_interval");

    // The frequency measurement limit is distinct from the input extent.
    const auto input_file = std::filesystem::absolute("ctint-contract-g0.h5").string();
    const auto write_input = [&](std::size_t count, bool nonfinite = false, unsigned flavors = 2) {
        alps::hdf5::archive archive(input_file, "w");
        std::vector<std::complex<double>> values(count, {0.0, -0.5});
        if (nonfinite) values[0] = {std::numeric_limits<double>::infinity(), 0.0};
        for (unsigned flavor = 0; flavor < flavors; ++flavor)
            archive["/G0_" + std::to_string(flavor)] << values;
    };
    auto file_run = configuration();
    file_run.input["atomic"] = false;
    file_run.input["g0"] = input_file;
    file_run.parameters["NMATSUBARA_MEASUREMENTS"] = 2;
    write_input(4);
    alps::ctint::prepare_run(file_run);
    require(file_run.parameters["EPSSQ_0"].as<double>() == 1. &&
            file_run.parameters["EPSSQ_1"].as<double>() == 1., "file-input moment defaults changed");
    write_input(3);
    rejects([&] { alps::ctint::prepare_run(file_run); }, "vector of length 4");
    {
        alps::hdf5::archive archive(input_file, "w");
        const std::vector<std::complex<double>> malformed(12, {0.0, 0.0});
        archive.write("/G0_0", malformed.data(), std::vector<std::size_t>{4, 3});
    }
    rejects([&] { alps::ctint::prepare_run(file_run); }, "complex vector");
    write_input(4, true);
    rejects([&] { alps::ctint::prepare_run(file_run); }, "must be finite");
    write_input(4);
    file_run.output["results"] = input_file;
    rejects([&] { alps::ctint::prepare_run(file_run); }, "replace input.g0");
    file_run.output["results"] = "ctint-contract-results.h5";
    file_run.parameters["FLAVORS"] = 4;
    rejects([&] { alps::ctint::prepare_run(file_run); }, "/G0_2");
    write_input(4, false, 4);
    alps::ctint::prepare_run(file_run);
    file_run.parameters["EPS_3"] = .25;
    alps::ctint::prepare_run(file_run);
    std::filesystem::remove(input_file);

    // Explicit thermalization must finish; prethermalization observations are excluded.
    InteractionExpansion simulation(run, 0);
    require(!simulation.is_thermalized(), "thermalized before any updates");
    simulation.update(); simulation.measure();
    require(simulation.collect_results().at("Sign").count() == 0, "warmup counted as measurement");
    simulation.update(); simulation.measure();
    require(simulation.is_thermalized() && simulation.fraction_completed() > 0.0,
            "thermalization never completed");
    simulation.update(); simulation.measure();
    require(simulation.fraction_completed() >= 1.0, "finite sweeps never complete");
    require(simulation.collect_results().at("Sign").count() == 2, "wrong measurement count");

    for (bool histogram : {false, true}) {
        run.parameters["HISTOGRAM_MEASUREMENT"] = histogram;
        alps::solvers::ctint(run);
        alps::hdf5::archive archive(run.output["results"].as<std::string>(), "r");
        alps::params parameters, execution;
        std::string application;
        archive["/parameters"] >> parameters;
        archive["/run_config/execution"] >> execution;
        archive["/run_config/application"] >> application;
        require(application == "ctint" && !parameters.exists("MAX_TIME"), "saved run sections");
        require(execution["seed"].as<int>() == 42, "saved execution provenance");
        alps::alea::hdf5_serializer codec(archive, "/simulation/results");
        auto observable = [&](std::string const& name, std::initializer_list<double> expected) {
            alps::alea::batch_result<double> result;
            deserialize(codec, archive.encode_segment(name), result);
            require(result.count() == 2 && result.size() == expected.size(), "atomic observable shape/count changed");
            std::size_t i = 0;
            for (auto value : expected)
                require(std::abs(result.mean()(i++) - value) < 1e-12,
                        "equal-time observable differs from atomic occupation-state oracle");
        };
        // Four equally probable occupation states (0, up, down, up+down).
        if (histogram) {
            observable("Sz_0", {0.});
            observable("Sz2_0", {.5});
            observable("Sz0_Sz0", {.5});
        } else {
            observable("densities", {.5, .5});
            observable("density_correlation", {.25});
            observable("n_i n_j", {.5, .25, .25, .5});
        }
        for (int flavor = 0; flavor < 2; ++flavor) {
            std::vector<std::complex<double>> frequency;
            std::vector<double> time;
            archive["/G_omega/" + std::to_string(flavor) + "/mean/value"] >> frequency;
            archive["/G_tau/" + std::to_string(flavor) + "/mean/value"] >> time;
            require(frequency.size() == 4 && time.size() == 9, "missing Green function output");
            for (std::size_t i = 0; i < frequency.size(); ++i)
                require(std::abs(frequency[i] - std::complex<double>(0., -2. / ((2. * i + 1.) * std::acos(-1.)))) < 1e-12,
                        "noninteracting frequency Green function differs from the analytic atom");
            for (auto value : time)
                require(std::abs(value + .5) < 1e-12,
                        "noninteracting time Green function differs from the analytic atom");
        }
    }
    std::filesystem::remove(run.output["results"].as<std::string>());
    require(!std::filesystem::exists("matrix_size"), "implicit text output created");
    require(!std::filesystem::exists("staggered_sz"), "implicit spin sidecar created");
}
