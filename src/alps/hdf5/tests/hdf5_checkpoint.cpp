// SPDX-License-Identifier: MIT
// Adapted from ALPS reconciliation 81e7ad00d, hdf5_lifecycle.cpp and
// mc/tests/checkpoint.cpp, for the current archive backend and GoogleTest.
#include <alps/hdf5/archive.hpp>
#include <alps/testing/temporary_directory.hpp>
#include <boost/filesystem.hpp>
#include <gtest/gtest.h>
#include <memory>
#include <stdexcept>

namespace {
struct HookFailure : std::runtime_error {
    HookFailure() : std::runtime_error("original checkpoint hook failure") {}
};
}

class CheckpointPublication : public ::testing::Test {
protected:
    alps::testing::TemporaryDirectory temporary;
    boost::filesystem::path filename = (temporary.path() / "state.h5").string();
    void seed() {
        alps::hdf5::archive ar(filename, "a");
        ar.write("/old", 17);
        ar.close();
    }
    void preserved() {
        alps::hdf5::archive ar(filename);
        int value = 0;
        ar.read("/old", value);
        EXPECT_EQ(value, 17);
        EXPECT_FALSE(ar.is_data("/new"));
    }
    void no_temporaries() {
        for (auto const & entry : boost::filesystem::directory_iterator(filename.parent_path()))
            EXPECT_EQ(entry.path().filename().string().find("state.h5.tmp."), std::string::npos);
    }
};

TEST_F(CheckpointPublication, SuccessReplacesTheCompleteFile) {
    seed();
    alps::hdf5::save_checkpoint(filename, [&](alps::hdf5::archive & ar) {
        preserved(); // The old file stays visible while serialization is running.
        ar.write("/new", 23);
    });
    alps::hdf5::archive ar(filename);
    int value = 0;
    ar.read("/new", value);
    EXPECT_EQ(value, 23);
    EXPECT_FALSE(ar.is_data("/old"));
    no_temporaries();
}

TEST_F(CheckpointPublication, HookFailurePreservesOldFileAndException) {
    seed();
    EXPECT_THROW(alps::hdf5::save_checkpoint(filename, [](alps::hdf5::archive & ar) {
        ar.write("/new", 23);
        throw HookFailure();
    }), HookFailure);
    preserved();
    no_temporaries();
}

TEST_F(CheckpointPublication, FailedFirstSavePublishesNothing) {
    EXPECT_THROW(alps::hdf5::save_checkpoint(filename, [](alps::hdf5::archive & ar) {
        ar.write("/new", 23);
        throw HookFailure();
    }), HookFailure);
    EXPECT_FALSE(boost::filesystem::exists(filename));
    no_temporaries();
}

TEST_F(CheckpointPublication, RetainedViewsCannotKeepAStagedFileOpen) {
    for (bool fail : {false, true}) {
        SCOPED_TRACE(fail);
        std::unique_ptr<alps::hdf5::archive> retained;
        auto save = [&](alps::hdf5::archive & ar) {
            ar.write("/new", 31);
            retained = std::make_unique<alps::hdf5::archive>(ar);
            if (fail) throw HookFailure();
        };
        if (fail) EXPECT_THROW(alps::hdf5::save_checkpoint(filename, save), HookFailure);
        else EXPECT_NO_THROW(alps::hdf5::save_checkpoint(filename, save));
        ASSERT_NE(retained, nullptr);
        EXPECT_FALSE(retained->is_open());
        int value = 0;
        EXPECT_THROW(retained->read("/new", value), alps::hdf5::archive_closed);
        EXPECT_THROW(retained->write("/late", 1), alps::hdf5::archive_closed);
        retained.reset();
        no_temporaries();
    }
}

TEST_F(CheckpointPublication, CallbackCanReleaseItsOwnView) {
    alps::hdf5::save_checkpoint(filename, [](alps::hdf5::archive & ar) {
        ar.write("/new", 37);
        ar.close();
    });
    alps::hdf5::archive ar(filename);
    int value = 0;
    ar.read("/new", value);
    EXPECT_EQ(value, 37);
    no_temporaries();
}

TEST_F(CheckpointPublication, RenameFailurePreservesDestinationAndCleansStage) {
    boost::filesystem::create_directory(filename);
    auto const marker = filename / "keep";
    boost::filesystem::create_directory(marker);
    EXPECT_THROW(alps::hdf5::save_checkpoint(filename, [](alps::hdf5::archive & ar) {
        ar.write("/new", 41);
    }), boost::filesystem::filesystem_error);
    EXPECT_TRUE(boost::filesystem::is_directory(marker));
    no_temporaries();
}
