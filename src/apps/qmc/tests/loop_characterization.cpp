// SPDX-License-Identifier: MIT
// Looper already exposes lattice-free temperature and sweep settings. Pin
// their distinct alias precedence, annealing, conversions and checkpoint rules
// alongside the scheduler-based solvers; do not unify their input contracts.
#include "../looper/looper/temperature.h"
#include "../looper/looper/montecarlo.h"
#include <alps/hdf5/archive.hpp>
#include <alps/osiris/xdrdump.h>
#include <alps/testing/temporary_directory.hpp>
#include <gtest/gtest.h>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <sstream>

namespace {
alps::Parameters parameters(const std::string& text) {
    alps::Parameters p;
    if (text.empty()) return p;
    std::istringstream input(text);
    p.parse(input);
    return p;
}

TEST(LoopTemperature, EvaluatesAliasesWithTemperatureFirst) {
    for (const auto& input : {"T=4*J;J=2", "TEMPERATURE=4*J;J=2",
                             "BETA=J/16;J=2", "INVERSE_TEMPERATURE=J/16;J=2"}) {
        SCOPED_TRACE(input);
        EXPECT_DOUBLE_EQ(looper::temperature(parameters(input))(), 8.);
    }
    EXPECT_DOUBLE_EQ(looper::temperature(parameters(
        "T=2;TEMPERATURE=3;BETA=4;INVERSE_TEMPERATURE=5;beta=6;Beta=7"))(), 2.);
    EXPECT_DOUBLE_EQ(looper::temperature(parameters(
        "TEMPERATURE=3;BETA=4;INVERSE_TEMPERATURE=5;beta=6"))(), 3.);
    EXPECT_DOUBLE_EQ(looper::temperature(parameters(
        "BETA=4;INVERSE_TEMPERATURE=5;beta=6"))(), .25);
}

TEST(LoopTemperature, IgnoresLowercaseAliasesAndDefersMissingTemperatureError) {
    for (const auto& input : {"", "beta=2", "Beta=2", "temperature=2"}) {
        SCOPED_TRACE(input);
        looper::temperature t(parameters(input));
        EXPECT_EQ(t.final(), -1.);
        EXPECT_THROW(t(), std::logic_error);
    }
    looper::temperature negative(parameters("T=-2"));
    EXPECT_EQ(negative.final(), -2.);
    EXPECT_THROW(negative(), std::logic_error);
    EXPECT_TRUE(std::isinf(looper::temperature(parameters("BETA=0"))()));
}

TEST(LoopTemperature, AnnealsAndStopsAtTheFirstIncompleteSegment) {
    looper::temperature t(parameters(
        "T=2;T_DURATION_0=4;T_START_0=10;T_DURATION_1=2;T_START_1=6;"
        "T_DURATION_2=3;T_DURATION_3=5;T_START_3=100"));
    EXPECT_EQ(t.initial(), 10.);
    EXPECT_EQ(t.final(), 2.);
    EXPECT_EQ(t.annealing_steps(), 6);
    EXPECT_THROW(t(), std::invalid_argument);
    EXPECT_DOUBLE_EQ(t(1), 10.);
    EXPECT_DOUBLE_EQ(t(3), 8.);
    EXPECT_DOUBLE_EQ(t(5), 6.);
    EXPECT_DOUBLE_EQ(t(6), 4.);
    EXPECT_DOUBLE_EQ(t(7), 2.);
    EXPECT_DOUBLE_EQ(t(20), 2.);
    EXPECT_THROW(t.set_beta(3.), std::invalid_argument);
}

TEST(LoopTemperature, ResetBetaUpdatesConstantSchedule) {
    looper::temperature t(parameters("T=2"));
    t.set_beta(4.);
    EXPECT_DOUBLE_EQ(t(), .25);
    // initial() retains the sequence entry, while operator() uses final_.
    EXPECT_DOUBLE_EQ(t.initial(), 2.);
}

TEST(LoopSweeps, DefaultsAndRangesUseOneEighthOfMinimum) {
    looper::mc_steps defaults(parameters(""));
    EXPECT_EQ(defaults.sweeps().min(), 65536u);
    EXPECT_EQ(defaults.sweeps().max(), std::numeric_limits<unsigned int>::max());
    EXPECT_EQ(defaults.thermalization(), 8192);
    looper::mc_steps ranged(parameters("SWEEPS=\"[17:23]\""));
    EXPECT_EQ(ranged.sweeps().min(), 17u);
    EXPECT_EQ(ranged.sweeps().max(), 23u);
    EXPECT_EQ(ranged.thermalization(), 2);
    for (int i = 0; i < 19; ++i) ++ranged;
    EXPECT_TRUE(ranged.is_thermalized());
    EXPECT_DOUBLE_EQ(ranged.progress(), 1.);
    EXPECT_TRUE(ranged.can_work());
    for (int i = 0; i < 6; ++i) ++ranged;
    EXPECT_FALSE(ranged.can_work());
}

TEST(LoopSweeps, ThermalizationEvaluatesExpressionsAndWrapsNegativeValues) {
    looper::mc_steps expression(parameters("SWEEPS=20;THERMALIZATION=2*J;J=3"));
    EXPECT_EQ(expression.thermalization(), 6);
    looper::mc_steps negative(parameters("SWEEPS=20;THERMALIZATION=-1"));
    EXPECT_EQ(negative.thermalization(), -1);
    for (int i = 0; i < 20; ++i) ++negative;
    EXPECT_FALSE(negative.is_thermalized());
    EXPECT_TRUE(negative.can_work());
    EXPECT_DOUBLE_EQ(negative.progress(), 20./19.);
}

TEST(LoopSweeps, CheckpointRestoresScheduleAndAllowsExplicitExtension) {
    looper::mc_steps original(parameters("SWEEPS=\"[17:23]\";THERMALIZATION=3"));
    for (int i = 0; i < 9; ++i) ++original;
    alps::testing::TemporaryDirectory directory;
    const auto file = (directory.path() / "loop.xdr").string();
    { alps::OXDRFileDump out(file); original.save(out); }
    looper::mc_steps resumed(parameters("SWEEPS=100;THERMALIZATION=10"));
    { alps::IXDRFileDump in(file); resumed.load(in); }
    EXPECT_EQ(resumed(), 9u);
    EXPECT_EQ(resumed.thermalization(), 3);
    EXPECT_EQ(resumed.sweeps().min(), 17u);
    EXPECT_EQ(resumed.sweeps().max(), 23u);
    EXPECT_DOUBLE_EQ(resumed.progress(), original.progress());
    resumed.set_sweeps(100);
    EXPECT_EQ(resumed.thermalization(), 3);
    EXPECT_DOUBLE_EQ(resumed.progress(), 9./103.);
    if (const char* dump = std::getenv("ALPS_QMC_DUMP_DIR")) {
        alps::OXDRFileDump out(std::string(dump) + "/loop_schedule.xdr");
        resumed.save(out);
        alps::hdf5::archive ar(std::string(dump) + "/loop_schedule.h5", "w");
        ar["/step"] << resumed();
        ar["/thermalization"] << resumed.thermalization();
        ar["/progress"] << resumed.progress();
    }
}
} // namespace
