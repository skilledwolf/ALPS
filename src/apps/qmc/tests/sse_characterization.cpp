// SPDX-License-Identifier: MIT
// Characterize the existing SSE scheduler settings before extracting them.
#include "../sse4/sse.h"
#include "construction_draws.h"

#include <alps/osiris/xdrdump.h>
#include <alps/testing/temporary_directory.hpp>
#include <gtest/gtest.h>

#include <cstdlib>
#include <cmath>
#include <filesystem>
#include <limits>
#include <string>

namespace {
alps::Parameters parameters() {
    alps::Parameters p;
    p["LATTICE"] = "chain lattice";
    p["MODEL"] = "spin";
    p["L"] = 4;
    p["local_S"] = 0.5;
    p["J"] = 1;
    p["T"] = 2;
    p["SWEEPS"] = 20;
    p["THERMALIZATION"] = 2;
    p["SEED"] = 713;
    p["NUMBER_OF_WORMS_PER_SWEEP"] = 2;
    return p;
}

alps::ProcessList processes() { return alps::ProcessList(1, alps::Process(0)); }

void steps(SSE_run& run, unsigned count) {
    for (unsigned i = 0; i < count; ++i)
        run.dostep();
}

// These optional artifacts are compared against a separately built baseline.
// The empty worker log avoids wall-clock fields. Assertions below deliberately
// cover settings only, so independent solver corrections need no new answers.
void dump(SSE_run& run, const std::string& name) {
    const char* root = std::getenv("ALPS_QMC_DUMP_DIR");
    if (!root || !*root)
        return;
    std::filesystem::create_directories(root);
    const auto stem = (std::filesystem::path(root) / ("sse_" + name)).string();
    {
        alps::hdf5::archive ar(stem + ".h5", "w");
        static_cast<const alps::scheduler::Worker&>(run).save(ar);
        ar["/results"] << run.measurements();
    }
    alps::OXDRFileDump checkpoint(stem + ".xdr");
    run.save(checkpoint);
}
} // namespace

TEST(SSECharacterization, TemperatureAliasesHaveExpressionPrecedence) {
    auto p = parameters();
    p["U"] = 0.25;
    const char* names[] = {"Beta", "beta", "BETA", "T", "TEMPERATURE", "temperature"};
    const char* values[] = {"2*U", "3*U", "4*U", "8*U", "16*U", "32*U"};
    const double expected[] = {0.5, 0.75, 1.0, 0.5, 0.25, 0.125};
    for (unsigned i = 0; i < 6; ++i)
        p[names[i]] = values[i];
    for (unsigned i = 0; i < 6; ++i) {
        SCOPED_TRACE(names[i]);
        SSE_run run(processes(), p, 0);
        EXPECT_DOUBLE_EQ(run.beta(), expected[i]);
        steps(run, 8);
        dump(run, "temperature_" + std::to_string(i));
        p.erase(names[i]);
    }
}

TEST(SSECharacterization, MissingAndNegativeTemperaturesAreRejected) {
    auto p = parameters();
    p.erase("T");
    EXPECT_THROW(SSE_run(processes(), p, 0), std::exception);
    p["BETA"] = -1;
    EXPECT_THROW(SSE_run(processes(), p, 0), std::out_of_range);
    p.erase("BETA");
    p["T"] = -2;
    EXPECT_THROW(SSE_run(processes(), p, 0), std::out_of_range);
}

TEST(SSECharacterization, SweepAliasesHaveFixedPrecedence) {
    auto p = parameters();
    p["THERMALIZATION"] = 0;
    p["SWEEPS"] = 8;
    p["MCS"] = 16;
    p["Steps"] = 32;
    const char* names[] = {"SWEEPS", "MCS", "Steps"};
    const double expected[] = {1.0 / 8, 1.0 / 16, 1.0 / 32};
    for (unsigned i = 0; i < 3; ++i) {
        SCOPED_TRACE(names[i]);
        SSE_run run(processes(), p, 0);
        EXPECT_TRUE(run.is_thermalized());
        run.dostep();
        EXPECT_DOUBLE_EQ(run.work_done(), expected[i]);
        steps(run, 7);
        dump(run, "sweeps_" + std::to_string(i));
        p.erase(names[i]);
    }
    EXPECT_THROW(SSE_run(processes(), p, 0), std::exception);
}

TEST(SSECharacterization, SweepCountsUsePlainUnsignedIntegerConversion) {
    auto p = parameters();
    p["SWEEPS"] = "2*10";
    EXPECT_THROW(SSE_run(processes(), p, 0), boost::bad_lexical_cast);
    p["SWEEPS"] = "20.5";
    EXPECT_THROW(SSE_run(processes(), p, 0), boost::bad_lexical_cast);
    p["SWEEPS"] = -1;
    p["THERMALIZATION"] = 0;
    SSE_run run(processes(), p, 0);
    run.dostep();
    EXPECT_DOUBLE_EQ(run.work_done(), 1.0 / double(std::numeric_limits<boost::uint64_t>::max()));
}

TEST(SSECharacterization, ZeroSweepTargetIsAcceptedAtConstruction) {
    auto p = parameters();
    p["SWEEPS"] = 0;
    p["THERMALIZATION"] = 0;
    SSE_run run(processes(), p, 0);
    EXPECT_TRUE(run.is_thermalized());
    EXPECT_TRUE(std::isnan(run.work_done()));
    run.dostep();
    EXPECT_TRUE(std::isinf(run.work_done()));
}

TEST(SSECharacterization, ThermalizationAliasesHaveFixedPrecedence) {
    auto p = parameters();
    p["THERMALIZATION"] = 3;
    p["thermalization"] = 1;
    SSE_run upper(processes(), p, 0);
    steps(upper, 2);
    EXPECT_FALSE(upper.is_thermalized());
    EXPECT_DOUBLE_EQ(upper.work_done(), 0.0);
    upper.dostep();
    EXPECT_TRUE(upper.is_thermalized());
    EXPECT_DOUBLE_EQ(upper.work_done(), 0.0);
    upper.dostep();
    EXPECT_DOUBLE_EQ(upper.work_done(), 1.0 / 20);
    dump(upper, "thermalization_upper");

    p.erase("THERMALIZATION");
    SSE_run lower(processes(), p, 0);
    EXPECT_FALSE(lower.is_thermalized());
    lower.dostep();
    EXPECT_TRUE(lower.is_thermalized());
    EXPECT_DOUBLE_EQ(lower.work_done(), 0.0);
    lower.dostep();
    EXPECT_DOUBLE_EQ(lower.work_done(), 1.0 / 20);
    dump(lower, "thermalization_lower");
}

TEST(SSECharacterization, DefaultThermalizationTruncatesOneTenthOfSweeps) {
    auto p = parameters();
    p.erase("THERMALIZATION");
    p["SWEEPS"] = 29;
    SSE_run run(processes(), p, 0);
    run.dostep();
    EXPECT_FALSE(run.is_thermalized());
    run.dostep();
    EXPECT_TRUE(run.is_thermalized());
    EXPECT_DOUBLE_EQ(run.work_done(), 0.0);
    run.dostep();
    EXPECT_DOUBLE_EQ(run.work_done(), 1.0 / 29);
    dump(run, "thermalization_default");
}

TEST(SSECharacterization, ThermalizationUsesPlainUnsignedIntegerConversion) {
    auto p = parameters();
    p["THERMALIZATION"] = "2*1";
    EXPECT_THROW(SSE_run(processes(), p, 0), boost::bad_lexical_cast);
    p["THERMALIZATION"] = "2.5";
    EXPECT_THROW(SSE_run(processes(), p, 0), boost::bad_lexical_cast);
    p["THERMALIZATION"] = -1;
    SSE_run run(processes(), p, 0);
    steps(run, 3);
    EXPECT_FALSE(run.is_thermalized());
    EXPECT_DOUBLE_EQ(run.work_done(), 0.0);
}

TEST(SSECharacterization, ThermalizationFractionUsesStrictComparison) {
    auto p = parameters();
    p["THERMALIZATION"] = 10;
    SSE_run run(processes(), p, 0);
    steps(run, 8);
    EXPECT_FALSE(run.is_thermalization_done(0.8));
    run.dostep();
    EXPECT_TRUE(run.is_thermalization_done(0.8));
    EXPECT_FALSE(run.is_thermalized());
}

TEST(SSECharacterization, SweepChangesAcceptPositiveAliasesAndRejectOtherSettings) {
    auto p = parameters();
    p["THERMALIZATION"] = 0;
    SSE_run run(processes(), p, 0);
    steps(run, 4);
    EXPECT_DOUBLE_EQ(run.work_done(), 4.0 / 20);
    for (const char* name : {"SWEEPS", "MCS", "Steps"}) {
        SCOPED_TRACE(name);
        EXPECT_TRUE(run.change_parameter(name, alps::StringValue(40)));
        EXPECT_DOUBLE_EQ(run.work_done(), 4.0 / 40);
        EXPECT_FALSE(run.change_parameter(name, alps::StringValue(0)));
        EXPECT_DOUBLE_EQ(run.work_done(), 4.0 / 40);
    }
    EXPECT_FALSE(run.change_parameter("THERMALIZATION", alps::StringValue(3)));
    EXPECT_FALSE(run.change_parameter("T", alps::StringValue(3)));
    EXPECT_FALSE(run.change_parameter("SEED", alps::StringValue(3)));
    // Worker::set_parameters separately allows SEED without reseeding.
    alps::Parameters updates;
    updates["SEED"] = 999;
    updates["SWEEPS"] = 80;
    EXPECT_NO_THROW(run.set_parameters(updates));
    EXPECT_DOUBLE_EQ(run.work_done(), 4.0 / 80);
    updates.clear();
    updates["THERMALIZATION"] = 3;
    EXPECT_THROW(run.set_parameters(updates), std::runtime_error);
    dump(run, "changed_sweeps");
}

TEST(SSECharacterization, ResumeRestoresThermalizationButUsesNewSweepCount) {
    auto p = parameters();
    p["THERMALIZATION"] = 3;
    SSE_run original(processes(), p, 0);
    steps(original, 5);
    EXPECT_DOUBLE_EQ(original.work_done(), 2.0 / 20);
    dump(original, "resume_before");
    alps::testing::TemporaryDirectory directory;
    const auto checkpoint = (directory.path() / "run.xdr").string();
    {
        alps::OXDRFileDump out(checkpoint);
        original.save(out);
    }
    p["THERMALIZATION"] = 9;
    p["SWEEPS"] = 40;
    SSE_run resumed(processes(), p, 0);
    EXPECT_FALSE(resumed.is_thermalized());
    {
        alps::IXDRFileDump in(checkpoint);
        resumed.load(in);
    }
    EXPECT_TRUE(resumed.is_thermalized());
    EXPECT_DOUBLE_EQ(resumed.work_done(), 2.0 / 40);
    resumed.dostep();
    EXPECT_DOUBLE_EQ(resumed.work_done(), 3.0 / 40);
    dump(resumed, "resume_after");
}

// The start configuration is drawn while the run is constructed: one number
// per site of the four-site chain.
TEST(SSECharacterization, DrawsTheStartConfigurationDuringConstruction) {
    SSE_run run(processes(), parameters(), 0);
    EXPECT_EQ(construction_draws(parameters(), [&] { return run.mrandom_real(); }), 4);
}

// Constructing a run reseeds the process-wide disorder generator from
// DISORDERSEED, so expressions evaluated beforehand, for example to check the
// parameters, do not change the disorder the run sees.
TEST(SSECharacterization, ConstructionReseedsTheDisorderGenerator) {
    alps::Parameters p = parameters();
    p["DISORDERSEED"] = 7;
    alps::Disorder::seed(7);
    const double first = alps::Disorder::random();
    const double a = alps::evaluate<double>("random()", p);
    const double b = alps::evaluate<double>("random()", p);
    EXPECT_NE(a, b);
    SSE_run run(processes(), p, 0);
    EXPECT_EQ(alps::Disorder::random(), first);
}
