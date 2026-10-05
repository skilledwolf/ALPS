// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#include <alps/mcbase.hpp>
#include <alps/mcmpiadapter.hpp>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
void require(bool condition, char const* message) {
    int success = condition;
    MPI_Allreduce(MPI_IN_PLACE, &success, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if (!success) throw std::runtime_error(message);
}
template<class F> void rejects(F operation, std::string const& message = {}) {
    bool failed = false;
    try { operation(); }
    catch (std::exception const& error) { failed = message.empty() || message == error.what(); }
    require(failed, "invalid collective MC operation was accepted or lost its error");
}
struct every_step {
    bool pending() const { return true; }
    void update(double) {}
};
class simulation : public alps::mcbase {
public:
    simulation(alps::params const& parameters, std::size_t rank)
        : mcbase(parameters, rank) {
        measurements.emplace("value", std::make_shared<alps::alea::batch_acc<double>>(1, rank ? 6 : 4, 2));
        measurements.emplace("empty", std::make_shared<alps::alea::batch_acc<double>>(2, 4));
    }
    void update() override {
        if (fail_update) throw std::runtime_error("rank update failure");
        ++steps;
    }
    void measure() override {
        if (fail_measure) throw std::runtime_error("rank measurement failure");
        *measurements.at("value") << alps::alea::make_adapter(double(steps));
    }
    double fraction_completed() const override { return completed ? 1. : steps/10.; }
    int steps = 0;
    bool fail_update = false, fail_measure = false, completed = false;
};
using mpi_simulation = alps::mcmpiadapter<simulation, every_step>;

void contract(boost::mpi::communicator const& comm) {
    auto const rank = comm.rank();
    for (int active : {0, 1}) {
        mpi_simulation sim({}, comm);
        if (rank == active)
            for (int i=1; i<=7; ++i)
                *sim.get_measurements().at("value") << alps::alea::make_adapter(double(i));
        auto const raw = sim.get_measurements().at("value")->result();
        auto result = sim.collect_results();
        require(sim.get_measurements().at("value")->result() == raw,
                "collective result collection changed the live accumulator");
        require(rank ? result.empty() : result.size() == 2, "results survived on wrong MPI rank");
        require(rank || (result.at("value").count() == 7 && result.at("value").count2() == 13.
            && result.at("value").mean()(0) == 4. && result.at("empty").count() == 0),
            "empty rank changed counts, partial-bin weights or means");
    }
    mpi_simulation sim({}, comm);
    for (int i=1; i<=(rank ? 5 : 7); ++i)
        *sim.get_measurements().at("value") << alps::alea::make_adapter(double(i+100*rank));
    auto result = sim.collect_results({"value"});
    bool correct = true;
    if (!rank) {
        double sums[]{3., 7., 11., 7., 203., 207., 105.};
        double counts[]{2., 2., 2., 1., 2., 2., 1.};
        auto const mean = 543./12.;
        double variance = 0.;
        for (int i=0; i<7; ++i) variance += counts[i]*std::pow(sums[i]/counts[i]-mean, 2);
        variance /= 12.-22./12.;
        auto const error = std::sqrt(variance*22./144.);
        auto const& pooled = result.at("value");
        correct = pooled.count() == 12 && pooled.count2() == 22.
            && pooled.mean()(0) == mean && std::abs(pooled.stderror()(0)-error) < 1e-12;
    }
    require(correct, "MPI pooling changed independent raw-bin estimates");
    rejects([&] { sim.collect_results(rank ? simulation::result_names_type{"empty"}
                                         : simulation::result_names_type{"value"}); });
    rejects([&] { sim.collect_results(rank ? simulation::result_names_type{"value", "empty"}
                                         : simulation::result_names_type{"value"}); });
    if (rank) sim.get_measurements().erase("value");
    rejects([&] { sim.collect_results({"value"}); });
    sim.get_measurements()["value"] = std::make_shared<alps::alea::batch_acc<double>>(rank ? 2 : 1, 4);
    rejects([&] { sim.collect_results({"value"}); });
    if (rank) sim.get_measurements()["value"].reset();
    rejects([&] { sim.collect_results({"value"}); });

    for (bool measurement : {false, true}) {
        mpi_simulation failing({}, comm);
        failing.fail_update = rank && !measurement;
        failing.fail_measure = rank && measurement;
        rejects([&] { failing.run([] { return false; }); },
                measurement ? "rank measurement failure" : "rank update failure");
    }
    mpi_simulation callback_failure({}, comm);
    rejects([&] { callback_failure.run([&] {
        if (rank) throw std::runtime_error("rank callback failure");
        return false;
    }); }, "rank callback failure");
    mpi_simulation stopped({}, comm);
    require(!stopped.run([&] { return rank != 0; }) && stopped.steps == 0,
            "initial MPI stop performed work or returned inconsistent completion");
    mpi_simulation completed({}, comm);
    completed.completed = true;
    require(completed.run([] { return false; }) && completed.steps == 0,
            "completed MPI simulation performed an extra update");
    mpi_simulation running({}, comm);
    require(running.run([] { return false; }) && running.steps == 5,
            "MPI aggregate progress changed simulation work allocation");
}
}

int main(int argc, char** argv) {
    boost::mpi::environment environment(argc, argv);
    boost::mpi::communicator communicator;
    try {
        require(communicator.size() == 2, "MC MPI contract requires two ranks");
        contract(communicator);
        if (!communicator.rank()) std::cout << "Native MC MPI contracts passed\n";
        return 0;
    } catch (std::exception const& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
