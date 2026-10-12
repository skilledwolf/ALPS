/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2012 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <gtest/gtest.h>
#include <alps/testing/temporary_directory.hpp>

#include <iostream>
#include <sstream>
#include <vector>
#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/vector.hpp>
#include <alps/ngs/cast.hpp>
#include <hdf5.h>

TEST(Hdf5, Valgrind) {
    alps::testing::TemporaryDirectory temporary;

    for (int i=0; i<100; ++i) {
        std::vector<double> vec(10,2.);
        alps::hdf5::archive ar((temporary.path() / "test_hdf5_valgrind.h5").string(), "w");
        std::ostringstream ss;
        ss << "/vec" << i;
        SCOPED_TRACE(i);
        ar << alps::make_pvp(ss.str(), vec);
        std::vector<double> restored;
        ar[ss.str()] >> restored;
        EXPECT_EQ(restored, vec);
    }

}

// Adapted from experiment fc02e70f0: keep the existing conversions, but test
// ownership on every completed VLEN read, including exception paths.
class StringReadOwnership : public ::testing::Test {
protected:
    alps::testing::TemporaryDirectory temporary;
    alps::hdf5::archive archive{(temporary.path() / "strings.h5").string(), "w"};
    std::string const overflow = std::string(4096, '9');
    std::vector<std::string> const values{"11", overflow, "33", "44", "55", "66"};
    std::vector<std::size_t> const shape{2, 3};

    void SetUp() override {
        archive.create_group("/group");
        archive.write("/scalar", overflow);
        archive.write("/group/@scalar", overflow);
        archive.write("/scalar/@scalar", overflow);
        archive.write("/vector", values.data(), shape);
        archive.write("/group/@vector", values.data(), shape);
        archive.write("/scalar/@vector", values.data(), shape);
    }

    void check_vector(char const* path, std::vector<std::size_t> const& chunk,
                      std::vector<std::size_t> const& offset,
                      std::vector<std::string> const& expected) {
        SCOPED_TRACE(path);
        const auto before = H5Fget_obj_count(H5F_OBJ_ALL, H5F_OBJ_ALL);
        ASSERT_GE(before, 0);
        std::vector<std::string> actual(expected.size());
        archive.read(path, actual.data(), chunk, offset);
        EXPECT_EQ(actual, expected);
        std::vector<int> converted(expected.size());
        // Legacy vector string-to-number reads reject char* conversion. Do
        // not change that policy while fixing the leaked VLEN allocations.
        try {
            archive.read(path, converted.data(), chunk, offset);
            FAIL() << "Expected a string-to-number conversion error";
        } catch (alps::bad_cast const& error) {
            EXPECT_EQ(std::string(error.what()).find("cannot cast from "), 0u);
        }
        EXPECT_EQ(H5Fget_obj_count(H5F_OBJ_ALL, H5F_OBJ_ALL), before);
    }
};

TEST_F(StringReadOwnership, DatasetStrings) {
    std::vector<std::string> const first{"11", overflow}, last{"55", "66"};
    for (int repeat = 0; repeat < 16; ++repeat) {
        check_vector("/vector", shape, {0, 0}, values);
        // File extent is six; the allocation to reclaim contains two values.
        check_vector("/vector", {1, 2}, {1, 1}, last);
        check_vector("/vector", {1, 2}, {0, 0}, first);
    }
}

TEST_F(StringReadOwnership, AttributeStrings) {
    for (int repeat = 0; repeat < 16; ++repeat) {
        for (char const* path : {"/group/@vector", "/scalar/@vector"}) {
            check_vector(path, shape, {0, 0}, values);
            const auto before = H5Fget_obj_count(H5F_OBJ_ALL, H5F_OBJ_ALL);
            std::vector<std::string> selected(2);
            EXPECT_THROW(archive.read(path, selected.data(), {1, 2}, {0, 0}), std::logic_error);
            EXPECT_EQ(H5Fget_obj_count(H5F_OBJ_ALL, H5F_OBJ_ALL), before);
        }
    }
}

TEST_F(StringReadOwnership, ScalarStrings) {
    for (char const* path : {"/scalar", "/group/@scalar", "/scalar/@scalar"}) {
        SCOPED_TRACE(path);
        std::string actual;
        archive.read(path, actual);
        EXPECT_EQ(actual, overflow);
        try {
            int value = 0;
            archive.read(path, value);
            FAIL() << "Expected an out-of-range scalar conversion";
        } catch (std::out_of_range const& error) {
            EXPECT_EQ(error.what(), "integer parameter out of range: " + overflow);
        }
    }
}
