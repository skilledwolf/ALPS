// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
// Physical paired-kink endpoint evidence: a stored0/1 must not undergo a
// second periodic normalization while a linked kink or checkpoint is copied.
#include "WKink.h"
#include <gtest/gtest.h>
#include <limits>
#include <vector>

namespace {
struct record_dump : alps::ODump {
    std::vector<double> times;
    void write_simple(int) override {}
    void write_simple(double value) override { times.push_back(value); }
    void write_simple(unsigned int) override {}
};
}
TEST(WormTimes, HistoricalNumericNormalizationRemainsUnchanged) {
    EXPECT_DOUBLE_EQ(double(time_struct(0.)), 1.);
    EXPECT_DOUBLE_EQ(double(time_struct(1.)), 0.);
    EXPECT_DOUBLE_EQ(double(time_struct(.25)), .25);
}
TEST(WormTimes, RawPhysicalEndpointsRemainEqualOnMovedAndNewLinkedKinks) {
    for (const double raw : {0., 1., .25}) {
        const auto time = time_struct::from_raw(raw);
        // One physical endpoint is moved; its partner is newly inserted.
        Kink<std::uint8_t> moved(.5, 0); moved.set_time(time); moved.set_id(17);
        auto partner = Kink<std::uint8_t>::at_time(time, 1); partner.set_id(17);
        EXPECT_EQ(moved.time(), partner.time());
        EXPECT_DOUBLE_EQ(double(moved.time()), raw);
        EXPECT_TRUE(moved.is_linked(partner));
        record_dump output; output << moved << partner;
        ASSERT_EQ(output.times.size(), 2); EXPECT_EQ(output.times[0], raw); EXPECT_EQ(output.times[1], raw);
        auto restored = Kink<std::uint8_t>::at_time(time_struct::from_raw(output.times[1]), partner.state());
        EXPECT_EQ(restored.time(), partner.time()); EXPECT_EQ(restored.state(), partner.state()); EXPECT_EQ(restored.id(), 0);
    }
    EXPECT_TRUE(std::signbit(double(time_struct::from_raw(-0.))));
}
TEST(WormTimes, InvalidRawTimesRejectInsteadOfPeriodicRepair) {
    for (const double raw : {-.1, 1.1, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()})
        EXPECT_THROW(time_struct::from_raw(raw), std::invalid_argument);
}
