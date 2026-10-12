// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
// Exercise the retained NGS adapter, independently of the modern ALEA reducer.
#include <alps/mcbase.hpp>
#include <alps/mcmpiadapter.hpp>
#include <gtest/gtest.h>

namespace {
class simulation : public alps::mcbase {
public:
    simulation(parameters_type const &p, std::size_t rank) : mcbase(p, rank) {}
    void update() override {}
    void measure() override {}
    double fraction_completed() const override { return 1.; }
};
using adapter = alps::mcmpiadapter<simulation>;

template<class Operation>
void rejected_on_all(const boost::mpi::communicator &world, Operation operation) {
    bool rejected = false;
    try { operation(); }
    catch (const alps::detail::collective_error &) { rejected = true; }
    EXPECT_TRUE(boost::mpi::all_reduce(world, rejected, std::logical_and<bool>()));
}

TEST(LegacyMPICollect, UnequalPositiveCountsKeepWeightedMean) {
    boost::mpi::communicator world;
    adapter sim({}, world);
    auto &obs = sim.get_measurements();
    obs << alps::accumulator::RealObservable("value")
        << alps::accumulator::RealVectorObservable("vector");
    const unsigned count = 16 * (world.rank() + 1);
    const double value = world.rank() + 1;
    for (unsigned i = 0; i < count; ++i) {
        obs["value"] << value;
        obs["vector"] << std::vector<double>{value, 2 * value};
    }
    auto results = sim.collect_results();
    if (world.rank() == 0) {
        // Counts are proportional to r; means are r. E[r] = sum r²/sum r.
        const double n = world.size();
        const double expected = (2 * n + 1) / 3;
        EXPECT_EQ(results["value"].count(), std::uint64_t(16 * n * (n + 1) / 2));
        EXPECT_NEAR(results["value"].mean<double>(), expected, 1e-12);
        EXPECT_NEAR(results["vector"].mean<std::vector<double>>()[0], expected, 1e-12);
        EXPECT_NEAR(results["vector"].mean<std::vector<double>>()[1], 2 * expected, 1e-12);
    }
}

TEST(LegacyMPICollect, AllEmptyRemainEmptyButMixedEmptyRejects) {
    boost::mpi::communicator world;
    adapter sim({}, world);
    sim.get_measurements() << alps::accumulator::RealObservable("value");
    auto results = sim.collect_results();
    EXPECT_EQ(results["value"].count(), 0u);
    if (world.rank() != world.size() - 1)
        for (unsigned i = 0; i < 16; ++i) sim.get_measurements()["value"] << 1.;
    rejected_on_all(world, [&] { sim.collect_results(); });
    // Agreement leaves the communicator usable; no reduction has started.
    EXPECT_EQ(boost::mpi::all_reduce(world, 1, std::plus<int>()), world.size());
}

TEST(LegacyMPICollect, NamesAndMissingObservablesRejectBeforeReduction) {
    boost::mpi::communicator world;
    adapter sim({}, world);
    sim.get_measurements() << alps::accumulator::RealObservable("value");
    rejected_on_all(world, [&] {
        sim.collect_results(world.rank() ? std::vector<std::string>{} : std::vector<std::string>{"value"});
    });
    rejected_on_all(world, [&] { sim.collect_results({"missing"}); });
}

TEST(LegacyMPICollect, ScalarVectorAndVectorExtentMismatchReject) {
    boost::mpi::communicator world;
    for (bool extent_only : {false, true}) {
        adapter sim({}, world);
        if (!extent_only && world.rank() == 0) {
            sim.get_measurements() << alps::accumulator::RealObservable("value");
            for (unsigned i = 0; i < 16; ++i) sim.get_measurements()["value"] << 1.;
        } else {
            sim.get_measurements() << alps::accumulator::RealVectorObservable("value");
            for (unsigned i = 0; i < 16; ++i)
                sim.get_measurements()["value"] << std::vector<double>(world.rank() ? 2 : 1, 1.);
        }
        rejected_on_all(world, [&] { sim.collect_results(); });
    }
}

TEST(LegacyMPICollect, BinConfigurationMismatchAndZeroReject) {
    boost::mpi::communicator world;
    for (bool zero : {false, true}) {
        alps::params p;
        p["BINNUMBER"] = zero ? 0 : (world.rank() ? 8 : 16);
        adapter sim(p, world);
        sim.get_measurements() << alps::accumulator::RealObservable("value");
        for (unsigned i = 0; i < 16; ++i) sim.get_measurements()["value"] << 1.;
        rejected_on_all(world, [&] { sim.collect_results(); });
    }
}
} // namespace

TEST(LegacyMPICollect, OptionalStatisticsMismatchRejects) {
    boost::mpi::communicator world;
    adapter sim({}, world);
    if (world.rank() == 0)
        sim.get_measurements() << alps::accumulator::RealVectorObservable("value");
    else
        sim.get_measurements() << alps::accumulator::SimpleRealVectorObservable("value");
    for (unsigned i = 0; i < 16; ++i)
        sim.get_measurements()["value"] << std::vector<double>{1., 2.};
    rejected_on_all(world, [&] { sim.collect_results(); });
}
