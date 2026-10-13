// SPDX-License-Identifier: MIT
// Unit tests for the task rules moved out of the scheduler: progress
// (progress.h), run checkpoint names (run_files.h) and result pooling
// (pooling.h). They use parameters and observables only; this file does not
// compile if one of the headers brings in a task, a worker or a scheduler.
#include <alps/scheduler/pooling.h>
#include <alps/scheduler/progress.h>
#include <alps/scheduler/run_files.h>

#if defined(ALPS_SCHEDULER_TASK_H) || defined(ALPS_SCHEDULER_WORKER_H) || \
    defined(ALPS_SCHEDULER_Scheduler_H) || defined(ALPS_SCHEDULER_MONTECARLO_H)
#error "the task rules must not depend on the scheduler"
#endif

#include <alps/alea.h>
#include <alps/parameter.h>
#include <gtest/gtest.h>

#include <sstream>
#include <string>
#include <vector>

namespace {

using alps::scheduler::measurement_pool;
using alps::scheduler::remaining_work;
using alps::scheduler::reported_percentage;
using alps::scheduler::run_checkpoint_name;
using alps::scheduler::work_factor;

alps::Parameters parameters(const std::string& text) {
    alps::Parameters p;
    std::istringstream input(text);
    p.parse(input);
    return p;
}

TEST(SchedulerTaskRules, WorkFactorIsAnExpressionOrOne) {
    EXPECT_DOUBLE_EQ(work_factor(parameters("SEED=1")), 1.);
    EXPECT_DOUBLE_EQ(work_factor(parameters("WORK_FACTOR=2.5")), 2.5);
    EXPECT_DOUBLE_EQ(work_factor(parameters("L=4;WORK_FACTOR=\"L*L\"")), 16.);
}

TEST(SchedulerTaskRules, RemainingWorkScalesTheWorkNotDone) {
    EXPECT_DOUBLE_EQ(remaining_work(parameters("SEED=1"), 0.), 1.);
    EXPECT_DOUBLE_EQ(remaining_work(parameters("SEED=1"), 0.25), 0.75);
    EXPECT_DOUBLE_EQ(remaining_work(parameters("WORK_FACTOR=4"), 0.25), 3.);
    EXPECT_DOUBLE_EQ(remaining_work(parameters("SEED=1"), 1.), 0.);
    EXPECT_DOUBLE_EQ(remaining_work(parameters("WORK_FACTOR=2"), 1.5), -1.);
}

TEST(SchedulerTaskRules, ReportedPercentageIsOneMinusWorkClipped) {
    EXPECT_DOUBLE_EQ(reported_percentage(0.75), 0.25);
    EXPECT_DOUBLE_EQ(reported_percentage(0.), 1.);
    EXPECT_DOUBLE_EQ(reported_percentage(-0.5), 1.);
    EXPECT_DOUBLE_EQ(reported_percentage(3.), 0.);
}

TEST(SchedulerTaskRules, NamesRunsAfterTheTaskFile) {
    EXPECT_EQ(run_checkpoint_name("job.task1.out.xml", {}), "job.task1.out.run1");
    EXPECT_EQ(run_checkpoint_name("/data/job.task1.xml", {"", "job.task1.run1"}),
              "job.task1.run2");
    EXPECT_EQ(run_checkpoint_name("job.task1.xml", {"job.task1.run2", "job.task1.run1"}),
              "job.task1.run3");
    EXPECT_EQ(run_checkpoint_name("job.task1.xml", {"job.task1.run2"}), "job.task1.run1");
    EXPECT_EQ(run_checkpoint_name("task", {}), "task.run1");
}

alps::ObservableSet run(double first, int count) {
    alps::ObservableSet set;
    alps::RealObservable value("Value");
    for (int i = 0; i < count; ++i)
        value << first + 0.1 * i;
    set << value;
    return set;
}

double mean(const alps::ObservableSet& set, const std::string& name) {
    return alps::RealObsevaluator(set[name]).mean();
}

double count(const alps::ObservableSet& set, const std::string& name) {
    return alps::RealObsevaluator(set[name]).count();
}

alps::ObservableSet merged(alps::ObservableSet a, const alps::ObservableSet& b) {
    a << b;
    return a;
}

TEST(SchedulerTaskRules, PoolsRunsIntoOneObservableInRunOrder) {
    const std::vector<alps::ObservableSet> r = {run(0., 2), run(4e15, 3), run(-3e15, 4),
                                                run(0., 5), run(1., 6)};
    measurement_pool pool;
    EXPECT_TRUE(pool.empty());
    for (const auto& set : r)
        pool.add(set);
    EXPECT_FALSE(pool.empty());
    const alps::ObservableSet pooled = pool.take(alps::ObservableSet());
    EXPECT_TRUE(pool.empty());

    const alps::ObservableSet expected =
        merged(merged(merged(r[0], r[1]), merged(r[2], r[3])), r[4]);
    EXPECT_EQ(pooled["Value"].number_of_runs(), 1u);
    EXPECT_EQ(count(pooled, "Value"), 20.);
    EXPECT_EQ(mean(pooled, "Value"), mean(expected, "Value"));

    measurement_pool reversed;
    for (auto it = r.rbegin(); it != r.rend(); ++it)
        reversed.add(*it);
    EXPECT_NE(mean(reversed.take(alps::ObservableSet()), "Value"), mean(pooled, "Value"));
}

TEST(SchedulerTaskRules, SimulationObservablesOnlyFillInMissingNames) {
    alps::ObservableSet simulation = run(100., 2);
    alps::RealObservable extra("Extra");
    extra << 1.;
    extra << 2.;
    simulation << extra;

    measurement_pool pool;
    pool.add(run(0., 3));
    const alps::ObservableSet pooled = pool.take(simulation);
    EXPECT_EQ(count(pooled, "Value"), 3.);
    EXPECT_DOUBLE_EQ(mean(pooled, "Value"), 0.1);
    EXPECT_DOUBLE_EQ(mean(pooled, "Extra"), 1.5);

    measurement_pool empty;
    EXPECT_TRUE(empty.take(simulation).empty());
}

} // namespace
