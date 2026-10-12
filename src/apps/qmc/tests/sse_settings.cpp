// SPDX-License-Identifier: MIT
#include "../sse4/settings.h"

#ifdef ALPS_SCHEDULER_WORKER_H
#error "Settings resolvers must not include scheduler workers"
#endif

#include <gtest/gtest.h>

#include <limits>

TEST(SSESettings, SweepAliasesHaveFixedPrecedence) {
    alps::Parameters p;
    p["SWEEPS"] = 8;
    p["MCS"] = 16;
    p["Steps"] = 32;
    EXPECT_EQ(dirloop_sse::resolve_sweeps(p), 8u);
    p.erase("SWEEPS");
    EXPECT_EQ(dirloop_sse::resolve_sweeps(p), 16u);
    p.erase("MCS");
    EXPECT_EQ(dirloop_sse::resolve_sweeps(p), 32u);
    p.erase("Steps");
    EXPECT_THROW(dirloop_sse::resolve_sweeps(p), std::exception);
}

TEST(SSESettings, UnselectedSweepAliasesAreNotConverted) {
    alps::Parameters p;
    p["SWEEPS"] = 8;
    p["MCS"] = "invalid";
    p["Steps"] = "invalid";
    EXPECT_EQ(dirloop_sse::resolve_sweeps(p), 8u);
    p.erase("SWEEPS");
    EXPECT_THROW(dirloop_sse::resolve_sweeps(p), boost::bad_lexical_cast);
    p["MCS"] = 16;
    EXPECT_EQ(dirloop_sse::resolve_sweeps(p), 16u);
}

TEST(SSESettings, SweepsPreserveUnsignedIntegerBoundaries) {
    alps::Parameters p;
    p["SWEEPS"] = 0;
    EXPECT_EQ(dirloop_sse::resolve_sweeps(p), 0u);
    p["SWEEPS"] = "18446744073709551615";
    EXPECT_EQ(dirloop_sse::resolve_sweeps(p), std::numeric_limits<boost::uint64_t>::max());
    p["SWEEPS"] = -1;
    EXPECT_EQ(dirloop_sse::resolve_sweeps(p), std::numeric_limits<boost::uint64_t>::max());
    p["SWEEPS"] = "18446744073709551616";
    EXPECT_THROW(dirloop_sse::resolve_sweeps(p), boost::bad_lexical_cast);
}

TEST(SSESettings, SweepsRequirePlainIntegers) {
    alps::Parameters p;
    p["U"] = 4;
    for (const char* value : {"2*U", "2.5", "", "invalid"}) {
        SCOPED_TRACE(value);
        p["SWEEPS"] = value;
        EXPECT_THROW(dirloop_sse::resolve_sweeps(p), boost::bad_lexical_cast);
    }
}

TEST(SSESettings, ThermalizationAliasesHaveFixedPrecedence) {
    alps::Parameters p;
    p["THERMALIZATION"] = 3;
    p["thermalization"] = 1;
    EXPECT_EQ(dirloop_sse::resolve_thermalization(p, 20), 3u);
    p.erase("THERMALIZATION");
    EXPECT_EQ(dirloop_sse::resolve_thermalization(p, 20), 1u);
    p.erase("thermalization");
    EXPECT_EQ(dirloop_sse::resolve_thermalization(p, 20), 2u);
}

TEST(SSESettings, UnselectedThermalizationAliasIsNotConverted) {
    alps::Parameters p;
    p["THERMALIZATION"] = 3;
    p["thermalization"] = "invalid";
    EXPECT_EQ(dirloop_sse::resolve_thermalization(p, 20), 3u);
    p.erase("THERMALIZATION");
    EXPECT_THROW(dirloop_sse::resolve_thermalization(p, 20), boost::bad_lexical_cast);
}

TEST(SSESettings, DefaultThermalizationUsesIntegerDivision) {
    alps::Parameters p;
    EXPECT_EQ(dirloop_sse::resolve_thermalization(p, 0), 0u);
    EXPECT_EQ(dirloop_sse::resolve_thermalization(p, 9), 0u);
    EXPECT_EQ(dirloop_sse::resolve_thermalization(p, 10), 1u);
    EXPECT_EQ(dirloop_sse::resolve_thermalization(p, 29), 2u);
    const auto max = std::numeric_limits<boost::uint64_t>::max();
    EXPECT_EQ(dirloop_sse::resolve_thermalization(p, max), max / 10);
}

TEST(SSESettings, ThermalizationPreservesUnsignedIntegerBoundaries) {
    alps::Parameters p;
    p["THERMALIZATION"] = 0;
    EXPECT_EQ(dirloop_sse::resolve_thermalization(p, 20), 0u);
    p["THERMALIZATION"] = "18446744073709551615";
    EXPECT_EQ(dirloop_sse::resolve_thermalization(p, 20), std::numeric_limits<boost::uint64_t>::max());
    p["THERMALIZATION"] = -1;
    EXPECT_EQ(dirloop_sse::resolve_thermalization(p, 20), std::numeric_limits<boost::uint64_t>::max());
    p["THERMALIZATION"] = "18446744073709551616";
    EXPECT_THROW(dirloop_sse::resolve_thermalization(p, 20), boost::bad_lexical_cast);
}

TEST(SSESettings, ThermalizationRequiresPlainIntegers) {
    alps::Parameters p;
    p["U"] = 4;
    for (const char* value : {"2*U", "2.5", "", "invalid"}) {
        SCOPED_TRACE(value);
        p["THERMALIZATION"] = value;
        EXPECT_THROW(dirloop_sse::resolve_thermalization(p, 20), boost::bad_lexical_cast);
    }
}
