// SPDX-License-Identifier: MIT
// Bespoke map regressions for the collection replacement contract adapted from
// ALPS reconciliation 81e7ad00d; retain this SDK's segment encoding.
#include <alps/hdf5/map.hpp>
#include <alps/testing/temporary_directory.hpp>
#include <gtest/gtest.h>
#include <map>
#include <string>

class Hdf5Map : public ::testing::Test {
protected:
    alps::testing::TemporaryDirectory temporary;
    alps::hdf5::archive ar{(temporary.path() / "maps.h5").string(), "a"};
};

TEST_F(Hdf5Map, ReplacesChildrenAndPersistsEmptyMaps) {
    std::map<std::string, int> source{{"old", 1}, {"removed", 2}};
    ar["/map"] << source;
    ar.write("/map/@owner", 17);
    ar.write("/sibling", 23);
    source = {{"A/B&C", 31}, {"..prefix", 37}};
    ar["/map"] << source;
    std::map<std::string, int> loaded{{"stale", 99}};
    ar["/map"] >> loaded;
    EXPECT_EQ(loaded, source);
    int value = 0;
    ar.read("/map/@owner", value); EXPECT_EQ(value, 17);
    source.clear();
    ar["/map"] << source;
    ASSERT_TRUE(ar.is_group("/map"));
    EXPECT_TRUE(ar.list_children("/map").empty());
    ar["/map"] >> loaded;
    EXPECT_TRUE(loaded.empty());
    ar.read("/map/@owner", value); EXPECT_EQ(value, 17);
    ar.read("/sibling", value); EXPECT_EQ(value, 23);
    ar.write("/scalar", 7);
    ar["/scalar"] << source;
    EXPECT_TRUE(ar.is_group("/scalar"));
}

TEST_F(Hdf5Map, InvalidKeysLeavePreviousContentsIntact) {
    ar["/map"] << std::map<std::string, int>{{"keep", 17}};
    for (auto const & key : std::vector<std::string>{"", ".", "..", "@attribute", "a@b", std::string("a\0b", 3)}) {
        SCOPED_TRACE(key);
        auto bad = std::map<std::string, int>{{key, 23}};
        EXPECT_THROW(ar["/map"] << bad, alps::hdf5::invalid_path);
        int value = 0;
        ar.read("/map/keep", value);
        EXPECT_EQ(value, 17);
    }
}

TEST_F(Hdf5Map, FailedLoadPreservesTheDestination) {
    ar.write("/map/first", 17);
    int const invalid[] = {1, 2};
    ar.write("/map/last", invalid, {2});
    std::map<std::string, int> destination{{"keep", 41}};
    auto const original = destination;
    EXPECT_THROW(ar["/map"] >> destination, std::exception);
    EXPECT_EQ(destination, original);
}

TEST_F(Hdf5Map, DuplicateDecodedKeysAreRejected) {
    ar.write("/map/&", 17);
    ar.write("/map/&#38;", 23);
    std::map<std::string, int> destination{{"keep", 41}};
    auto const original = destination;
    EXPECT_THROW(ar["/map"] >> destination, alps::hdf5::invalid_path);
    EXPECT_EQ(destination, original);
}
