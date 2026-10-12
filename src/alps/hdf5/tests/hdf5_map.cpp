// SPDX-License-Identifier: MIT
// Bespoke map regressions for the collection replacement contract adapted from
// ALPS reconciliation 81e7ad00d; retain this SDK's segment encoding.
#include <alps/hdf5/map.hpp>
#include <alps/testing/temporary_directory.hpp>
#include <gtest/gtest.h>
#include <hdf5.h>
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

TEST(Hdf5MapNativeChildren, RejectedBeforeReplacementOrLoading) {
    alps::testing::TemporaryDirectory temporary;
    for (std::string const kind : {"dangling", "external", "type", "zz@attribute", ".."}) {
        SCOPED_TRACE(kind);
        auto const filename = (temporary.path() / (kind + ".h5")).string();
        const std::string name = kind == "zz@attribute" || kind == ".." ? kind : "zz_unsupported";
        auto const child = "/map/" + name;
        {
            alps::hdf5::archive ar(filename, "a");
            ar["/map"] << std::map<std::string, int>{{"keep", 17}};
            ar.write("/map/@owner", 23);
            ar.write("/sibling", 31);
        }
        auto file = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
        ASSERT_GE(file, 0);
        if (kind == "dangling") {
            ASSERT_GE(H5Lcreate_soft("/missing", file, child.c_str(), H5P_DEFAULT, H5P_DEFAULT), 0);
        } else if (kind == "external") {
            ASSERT_GE(H5Lcreate_external("missing.h5", "/value", file, child.c_str(),
                                        H5P_DEFAULT, H5P_DEFAULT), 0);
        } else if (kind == "type") {
            auto type = H5Tcopy(H5T_NATIVE_INT);
            ASSERT_GE(type, 0);
            ASSERT_GE(H5Tcommit2(file, child.c_str(), type, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), 0);
            ASSERT_GE(H5Tclose(type), 0);
        } else {
            auto space = H5Screate(H5S_SCALAR);
            auto data = H5Dcreate2(file, child.c_str(), H5T_NATIVE_INT, space,
                                  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            ASSERT_GE(data, 0);
            const int value = 43;
            ASSERT_GE(H5Dwrite(data, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value), 0);
            ASSERT_GE(H5Dclose(data), 0);
            ASSERT_GE(H5Sclose(space), 0);
        }
        ASSERT_GE(H5Fclose(file), 0);
        {
            alps::hdf5::archive ar(filename, "a");
            const std::map<std::string, int> replacement{{"new", 37}};
            EXPECT_THROW(ar["/map"] << replacement, alps::hdf5::invalid_path);
            int value = 0;
            ar.read("/map/keep", value); EXPECT_EQ(value, 17);
            ar.read("/map/@owner", value); EXPECT_EQ(value, 23);
            ar.read("/sibling", value); EXPECT_EQ(value, 31);
            EXPECT_FALSE(ar.is_data("/map/new"));
            auto destination = replacement;
            EXPECT_THROW(ar["/map"] >> destination, alps::hdf5::invalid_path);
            EXPECT_EQ(destination, replacement);
        }
        file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        ASSERT_GE(file, 0);
        EXPECT_GT(H5Lexists(file, child.c_str(), H5P_DEFAULT), 0);
        ASSERT_GE(H5Fclose(file), 0);
    }
}

TEST_F(Hdf5Map, RejectsDatasetAncestorsWithoutReplacingThem) {
    ar.write("/ancestor", 17);
    ar.write("/ancestor/@owner", 23);
    const std::map<std::string, int> replacement{{"new", 37}};
    EXPECT_THROW(ar["/ancestor/nested/map"] << replacement, alps::hdf5::invalid_path);
    int value = 0;
    ar.read("/ancestor", value); EXPECT_EQ(value, 17);
    ar.read("/ancestor/@owner", value); EXPECT_EQ(value, 23);
    EXPECT_FALSE(ar.is_group("/ancestor"));
}

TEST(Hdf5MapNativeTargets, RejectsDanglingTargetAndAncestorWithoutMutation) {
    alps::testing::TemporaryDirectory temporary;
    const auto filename = (temporary.path() / "targets.h5").string();
    {
        alps::hdf5::archive ar(filename, "a");
        ar.write("/sibling", 17);
    }
    auto file = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    ASSERT_GE(file, 0);
    ASSERT_GE(H5Lcreate_soft("/missing", file, "/dangling", H5P_DEFAULT, H5P_DEFAULT), 0);
    ASSERT_GE(H5Fclose(file), 0);
    {
        alps::hdf5::archive ar(filename, "a");
        const std::map<std::string, int> replacement{{"new", 37}};
        EXPECT_THROW(ar["/dangling"] << replacement, alps::hdf5::invalid_path);
        EXPECT_THROW(ar["/dangling/nested/map"] << replacement, alps::hdf5::invalid_path);
        int value = 0;
        ar.read("/sibling", value); EXPECT_EQ(value, 17);
        EXPECT_FALSE(ar.is_group("/missing"));
    }
    file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    ASSERT_GE(file, 0);
    EXPECT_GT(H5Lexists(file, "/dangling", H5P_DEFAULT), 0);
    ASSERT_GE(H5Fclose(file), 0);
}
