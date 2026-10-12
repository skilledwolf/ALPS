// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
// Independent evidence that retained worm kink/link state is initialized
// before copying, incrementing or legacy serialization; original dump layout
// and loaded identifiers remain unchanged.
#include "WRun.h"
#include <gtest/gtest.h>
#include <array>
#include <cstring>
#include <new>

namespace {
struct record_dump : alps::ODump {
    std::vector<unsigned int> unsigned_words;
    void write_simple(int) override {}
    void write_simple(double) override {}
    void write_simple(unsigned int value) override { unsigned_words.push_back(value); }
};
alps::Parameters configuration() {
    alps::Parameters input;
    input["MODEL"] = "spin"; input["LATTICE"] = "dimer";
    input["LATTICE_LIBRARY"] = std::string(ALPS_WORM_LATTICES);
    input["MODEL_LIBRARY"] = std::string(ALPS_WORM_MODELS);
    input["T"] = 1.2; input["SWEEPS"] = 400; input["THERMALIZATION"] = 0;
    input["SEED"] = 173; input["RNG"] = "mt19937"; input["J"] = 1;
    input["DISORDERSEED"] = 173;
    return input;
}
struct patterned_run {
    alignas(WRun) std::array<unsigned char, sizeof(WRun)> storage;
    WRun *run;
    patterned_run(const alps::Parameters &input, unsigned char pattern) {
        std::memset(storage.data(), pattern, storage.size());
        run = new (storage.data()) WRun(alps::ProcessList(1), input, 0);
    }
    ~patterned_run() { run->~WRun(); }
};
TEST(WormInitialState, KinkPayloadAndUnlinkedIdentityAreDefinedBeforeCopying) {
    Kink<std::uint8_t> empty, explicit_state(.25, 1);
    EXPECT_EQ(empty.state(), 0); EXPECT_EQ(empty.id(), 0);
    EXPECT_EQ(explicit_state.state(), 1); EXPECT_DOUBLE_EQ(explicit_state.time(), .25);
    EXPECT_EQ(explicit_state.id(), 0);
    auto copy = explicit_state;
    copy.set_id(13);
    EXPECT_EQ(copy.id(), 13); EXPECT_EQ(explicit_state.id(), 0);
    EXPECT_EQ(copy.state(), 1); EXPECT_DOUBLE_EQ(copy.time(), .25);
}
TEST(WormInitialState, FreshLegacyDumpCounterIsZeroIndependentOfPriorObjectStorage) {
    const auto input = configuration();
    for (const auto pattern : {0x55, 0xaa}) {
        patterned_run instance(input, static_cast<unsigned char>(pattern));
        record_dump output; instance.run->save(output);
        ASSERT_FALSE(output.unsigned_words.empty());
        // Existing v100 WRun::save writes last_id_ as its final unsigned field.
        EXPECT_EQ(output.unsigned_words.back(), 0);
    }
}
TEST(WormInitialState, FixedSeedActualEnergyTraceIsIndependentOfPriorObjectStorage) {
    const auto input = configuration();
    patterned_run first(input, 0x55), second(input, 0xaa);
    first.run->start_worker(); second.run->start_worker();
    for (unsigned step = 0; step < 100; ++step) {
        first.run->run(); second.run->run();
        const auto &a = first.run->get_measurements().get<alps::RealObservable>("Energy");
        const auto &b = second.run->get_measurements().get<alps::RealObservable>("Energy");
        EXPECT_EQ(a.count(), step + 1); EXPECT_EQ(b.count(), step + 1);
        EXPECT_DOUBLE_EQ(a.mean(), b.mean());
        record_dump first_state, second_state;
        first.run->save(first_state); second.run->save(second_state);
        EXPECT_EQ(first_state.unsigned_words, second_state.unsigned_words);
    }
}
} // namespace
