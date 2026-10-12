// SPDX-License-Identifier: MIT
// Adapted from ALPS reconciliation 81e7ad00d, mc/tests/checkpoint.cpp.
// Preserve the current serializers and test exact native RNG continuation.
#include <alps/mcbase.hpp>
#include <alps/testing/temporary_directory.hpp>
#include <boost/filesystem.hpp>
#include <gtest/gtest.h>
#include <stdexcept>

namespace {
struct HookFailure : std::runtime_error {
    HookFailure() : std::runtime_error("original simulation save failure") {}
};
class Simulation : public alps::mcbase {
public:
    Simulation() : mcbase(alps::params{}) {
        get_parameters()["SEED"] = 42;
        get_measurements() << alps::accumulator::RealObservable("value");
        get_measurements()["value"] << 1.;
    }
    using mcbase::save;
    using mcbase::load;
    void update() override {}
    void measure() override {}
    double fraction_completed() const override { return 1.; }
    double draw() { return random(); }
    void save(alps::hdf5::archive & ar) const override {
        mcbase::save(ar);
        ar["state"] << state;
        if (fail) throw HookFailure();
    }
    void load(alps::hdf5::archive & ar) override {
        mcbase::load(ar);
        ar["state"] >> state;
    }
    int state = 0;
    bool fail = false;
};
}

class McCheckpoint : public ::testing::Test {
protected:
    alps::testing::TemporaryDirectory temporary;
    boost::filesystem::path filename = (temporary.path() / "simulation.h5").string();
    void no_temporaries() {
        for (auto const & entry : boost::filesystem::directory_iterator(filename.parent_path()))
            EXPECT_EQ(entry.path().filename().string().find("simulation.h5.tmp."), std::string::npos);
    }
};

TEST_F(McCheckpoint, CompleteReplacementAndExactRandomContinuation) {
    { alps::hdf5::archive ar(filename, "a"); ar.write("/stale", 1); }
    Simulation simulation;
    for (int i = 0; i < 5; ++i) simulation.draw();
    simulation.state = 17;
    simulation.save(filename);
    Simulation restored;
    restored.load(filename);
    EXPECT_EQ(restored.state, 17);
    for (int i = 0; i < 16; ++i) EXPECT_EQ(restored.draw(), simulation.draw());
    alps::hdf5::archive ar(filename);
    EXPECT_FALSE(ar.is_data("/stale"));
    EXPECT_TRUE(ar.is_data("/simulation/realizations/0/clones/0/state"));
    EXPECT_TRUE(ar.is_group("/parameters"));
    no_temporaries();
}

TEST_F(McCheckpoint, FailedSavePreservesLastCheckpoint) {
    Simulation simulation;
    simulation.state = 17;
    simulation.save(filename);
    simulation.state = 23;
    simulation.fail = true;
    EXPECT_THROW(simulation.save(filename), HookFailure);
    Simulation restored;
    restored.load(filename);
    EXPECT_EQ(restored.state, 17);
    no_temporaries();
}

TEST_F(McCheckpoint, FailedFirstSavePublishesNothing) {
    Simulation simulation;
    simulation.fail = true;
    EXPECT_THROW(simulation.save(filename), HookFailure);
    EXPECT_FALSE(boost::filesystem::exists(filename));
    no_temporaries();
}

TEST_F(McCheckpoint, PublicationFailurePreservesDestination) {
    boost::filesystem::create_directory(filename);
    auto const marker = filename / "keep";
    boost::filesystem::create_directory(marker);
    Simulation simulation;
    EXPECT_THROW(simulation.save(filename), boost::filesystem::filesystem_error);
    EXPECT_TRUE(boost::filesystem::is_directory(marker));
    no_temporaries();
}
