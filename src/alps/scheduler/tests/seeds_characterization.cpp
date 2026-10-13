// SPDX-License-Identifier: MIT
// Characterize how the legacy scheduler seeds random number generators: the
// engine of each worker, the process-wide disorder generator used by
// parameter expressions, and the SEED a task hands to each new run. These
// values decide which random numbers a run sees, so a refactor of the seed
// handling must keep them, including the conversions between signed and
// unsigned integers.
#include <alps/expression.h>
#include <alps/osiris/comm.h>
#include <alps/parameter.h>
#include <alps/scheduler/montecarlo.h>
#include <alps/scheduler/options.h>
#include <alps/scheduler/scheduler.h>
#include <gtest/gtest.h>

#include <sstream>
#include <stdexcept>
#include <string>

namespace {

alps::Parameters parameters(const std::string& text) {
    alps::Parameters p;
    std::istringstream input(text);
    p.parse(input);
    return p;
}

class probe_run : public alps::scheduler::MCRun {
public:
    probe_run(const alps::ProcessList& where, const alps::Parameters& p, int node)
        : alps::scheduler::MCRun(where, p, node) {}
    void dostep() {}
    double work_done() const { return 0.; }
    double draw() { return random_01(); }
    const alps::Parameters& parameters() const { return parms; }
};

const alps::ProcessList here(1, alps::local_process());

double first_draw(const std::string& text) {
    probe_run run(here, parameters(text), 0);
    return run.draw();
}

std::string construction_error(const std::string& text) {
    try {
        probe_run run(here, parameters(text), 0);
    } catch (const std::exception& error) {
        return error.what();
    }
    return "";
}

double disorder_draw_after_seed(unsigned int seed) {
    alps::Disorder::seed(seed);
    return alps::Disorder::random();
}

TEST(SchedulerSeeds, SeedsTheWorkerEngineWithSeedAsUnsignedInteger) {
    EXPECT_EQ(first_draw("SEED=5"), first_draw("SEED=5"));
    EXPECT_NE(first_draw("SEED=5"), first_draw("SEED=6"));
    EXPECT_EQ(first_draw("SEED=-5"), first_draw("SEED=4294967291"));
    EXPECT_EQ(construction_error("SEED=4294967296"),
              construction_error("T=1"));
    EXPECT_NE(construction_error("T=1"), "");
}

TEST(SchedulerSeeds, SeedsTheDisorderGeneratorWithDisorderseed) {
    const double after_9 = disorder_draw_after_seed(9);
    const double after_0 = disorder_draw_after_seed(0);
    const double after_max = disorder_draw_after_seed(4294967295u);

    alps::Disorder::seed(123);
    probe_run with_seed(here, parameters("SEED=1;DISORDERSEED=9"), 0);
    EXPECT_EQ(alps::Disorder::random(), after_9);

    alps::Disorder::seed(123);
    probe_run without_seed(here, parameters("SEED=1"), 0);
    EXPECT_EQ(alps::Disorder::random(), after_0);

    alps::Disorder::seed(123);
    probe_run negative_seed(here, parameters("SEED=1;DISORDERSEED=-1"), 0);
    EXPECT_EQ(alps::Disorder::random(), after_max);
}

// Expressions reseed the disorder generator from DISORDERSEED read as the
// worker reads it, so every seed the worker accepts works in expressions.
TEST(SchedulerSeeds, ExpressionsAcceptTheDisorderSeedsTheWorkerAccepts) {
    const auto p = parameters("SEED=1;DISORDERSEED=3000000000;J=2");
    probe_run run(here, p, 0);
    EXPECT_DOUBLE_EQ(alps::evaluate<double>("J", p), 2.);
    EXPECT_DOUBLE_EQ(alps::evaluate<double>("J", parameters("DISORDERSEED=7;J=2")), 2.);

    // A DISORDERSEED other than the last seed reseeds the generator.
    const double after_large = disorder_draw_after_seed(3000000001u);
    alps::Disorder::seed(123);
    alps::evaluate<double>("J", parameters("DISORDERSEED=3000000001;J=2"));
    EXPECT_EQ(alps::Disorder::random(), after_large);

    const double after_max = disorder_draw_after_seed(4294967295u);
    alps::Disorder::seed(123);
    alps::evaluate<double>("J", parameters("DISORDERSEED=-1;J=2"));
    EXPECT_EQ(alps::Disorder::random(), after_max);

    // The last seed again does not reseed: the stream continues.
    alps::Disorder::seed(3000000001u);
    alps::Disorder::random();
    const double second = alps::Disorder::random();
    alps::Disorder::seed(3000000001u);
    alps::Disorder::random();
    alps::evaluate<double>("J", parameters("DISORDERSEED=3000000001;J=2"));
    EXPECT_EQ(alps::Disorder::random(), second);
}

class TaskSeeds : public ::testing::Test {
protected:
    TaskSeeds() : scheduler(options, factory) {}

    // The SEED a single-process task keeps after creating its run, and the
    // SEED its run was created with.
    std::pair<std::string, std::string> seeds(const std::string& text) {
        alps::scheduler::MCSimulation task(here, parameters(text)); // constructs its runs
        const auto& run = dynamic_cast<const probe_run&>(*task.runs.at(0));
        return {std::string(task.get_parameters()["SEED"]),
                std::string(run.parameters()["SEED"])};
    }

    std::string task_error(const std::string& text) {
        try {
            seeds(text);
        } catch (const std::exception& error) {
            return error.what();
        }
        return "";
    }

    alps::scheduler::SimpleMCFactory<probe_run> factory;
    alps::scheduler::NoJobfileOptions options;
    alps::scheduler::Scheduler scheduler;
};

// Each run created advances the task's SEED by the processes per run. The
// sum is computed as int32_t(SEED) + uint32_t(cpus), in unsigned arithmetic.
TEST_F(TaskSeeds, AdvancesTheTaskSeedForEachRun) {
    EXPECT_EQ(seeds("SEED=5"), std::make_pair(std::string("6"), std::string("5")));
    EXPECT_EQ(seeds("T=1"), std::make_pair(std::string("1"), std::string("0")));
    EXPECT_EQ(seeds("SEED=2147483647"),
              std::make_pair(std::string("2147483648"), std::string("2147483647")));
    EXPECT_EQ(seeds("SEED=-5"), std::make_pair(std::string("4294967292"), std::string("-5")));
    EXPECT_EQ(seeds("SEED=7;DISORDERSEED=3").first, "8");
}

// A SEED above INT32_MAX seeds the first run's engine, but cannot be advanced.
TEST_F(TaskSeeds, CannotAdvanceSeedsAboveInt32Max) {
    EXPECT_NE(task_error("SEED=2147483648"), "");
    EXPECT_EQ(task_error("SEED=2147483648"), task_error("SEED=4294967295"));
    EXPECT_EQ(construction_error("SEED=2147483648"), "");
}

} // namespace
