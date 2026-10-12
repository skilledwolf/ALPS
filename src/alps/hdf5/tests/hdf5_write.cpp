// SPDX-License-Identifier: MIT
// Adapted from ALPS reconciliation 649d6c5a9, hdf5_write.cpp.
// Split into isolated GoogleTest cases for the current native test runner.
#include <gtest/gtest.h>
#include <alps/testing/temporary_directory.hpp>
#include <alps/hdf5/archive.hpp>

#include <hdf5.h>

#include <array>
#include <cstdio>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
void require(bool condition) {
    EXPECT_TRUE(condition);
}

template <typename Write>
void rejects(Write write) {
    auto const objects_before = H5Fget_obj_count(H5F_OBJ_ALL, H5F_OBJ_ALL);
    bool rejected = false;
    try {
        write();
    } catch (std::exception const &) {
        rejected = true;
    }
    require(rejected);
    require(H5Fget_obj_count(H5F_OBJ_ALL, H5F_OBJ_ALL) == objects_before);
}

template <typename T>
void array_writes(alps::hdf5::archive & ar, std::string const & name,
                  std::array<T, 6> const & values, std::array<T, 2> const & patch) {
    std::string const dataset = "/" + name;
    ar.write(dataset, values.data(), {2, 3});
    require(ar.dimensions(dataset) == 2);
    require(ar.extent(dataset) == std::vector<std::size_t>({2, 3}));
    std::array<T, 6> actual;
    ar.read(dataset, actual.data(), {2, 3});
    require(actual == values);

    ar.write(dataset + "/@kept", 42);
    ar.write(dataset, patch.data(), {2, 3}, {2, 1}, {0, 1});
    auto expected = values;
    expected[1] = patch[0];
    expected[4] = patch[1];
    ar.read(dataset, actual.data(), {2, 3});
    require(actual == expected);
    int kept = 0;
    ar.read(dataset + "/@kept", kept);
    require(kept == 42);
    ar.write(dataset, values.data(), {2, 3});
    ar.read(dataset, actual.data(), {2, 3});
    require(actual == values);
    ar.read(dataset + "/@kept", kept);
    require(kept == 42);

    // Attributes have no hyperslab transfer: complete replacement is valid,
    // whereas an attempted partial update must preserve the old value.
    for (auto const & attribute : {"/group/@" + name, dataset + "/@array"}) {
        ar.write(attribute, values.data(), {2, 3});
        require(ar.extent(attribute) == std::vector<std::size_t>({2, 3}));
        ar.read(attribute, actual.data(), {2, 3});
        require(actual == values);
        rejects([&] { ar.write(attribute, patch.data(), {2, 3}, {2, 1}, {0, 1}); });
        ar.read(attribute, actual.data(), {2, 3});
        require(actual == values);
    }
}

void attribute_replacements(alps::hdf5::archive & ar, std::string const & path) {
    std::array<int, 6> const values{1, 2, 3, 4, 5, 6};
    ar.write(path, 7);
    require(ar.is_scalar(path));
    ar.write(path, values.data(), {2, 3});
    require(!ar.is_scalar(path));
    require(ar.extent(path) == std::vector<std::size_t>({2, 3}));
    std::array<int, 6> actual;
    ar.read(path, actual.data(), {2, 3});
    require(actual == values);
    ar.write(path, values.data(), {3, 2});
    require(ar.extent(path) == std::vector<std::size_t>({3, 2}));
    ar.read(path, actual.data(), {3, 2});
    require(actual == values);
    ar.write(path, std::string("new type"));
    std::string text;
    ar.read(path, text);
    require(ar.is_scalar(path) && text == "new type");
    ar.write(path, 9);
    int scalar = 0;
    ar.read(path, scalar);
    require(ar.is_scalar(path) && scalar == 9);
    ar.write(path, static_cast<int const *>(nullptr), {});
    require(ar.is_null(path));
    ar.write(path, values.data(), {2, 3});
    ar.read(path, actual.data(), {2, 3});
    require(actual == values);
}

void replacements(alps::hdf5::archive & ar) {
    std::array<int, 6> const values{1, 2, 3, 4, 5, 6};
    ar.write("/replace", 1);
    ar.write("/replace/@kept", 42);
    ar.write("/replace", 2);
    int scalar = 0;
    ar.read("/replace/@kept", scalar);
    require(scalar == 42);
    ar.write("/replace", values.data(), {2, 3});
    require(!ar.is_attribute("/replace/@kept"));
    ar.write("/replace/@shape", 43);
    ar.write("/replace", values.data(), {3, 2});
    require(!ar.is_attribute("/replace/@shape"));
    ar.write("/replace/@class", 44);
    ar.write("/replace", 3);
    require(!ar.is_attribute("/replace/@class"));
    ar.write("/replace/@type", 45);
    ar.write("/replace", 4.5);
    require(!ar.is_attribute("/replace/@type"));
    double number = 0;
    ar.read("/replace", number);
    require(number == 4.5);

    ar.write("/group_replace/child", 91);
    ar.write("/group_replace/@old", 92);
    ar.write("/group_replace", values.data(), {2, 3});
    require(ar.is_data("/group_replace"));
    require(!ar.is_attribute("/group_replace/@old"));
    std::array<int, 6> actual;
    ar.read("/group_replace", actual.data(), {2, 3});
    require(actual == values);
}

void empty_writes(alps::hdf5::archive & ar) {
    for (auto path : {"/null", "/group/@null"}) {
        ar.write(path, static_cast<int const *>(nullptr), {});
        require(ar.is_null(path));
    }
    // A ranked empty array has a SIMPLE dataspace and retains its shape.
    for (auto path : {"/empty", "/group/@empty"}) {
        ar.write(path, static_cast<int const *>(nullptr), {2, 0});
        require(!ar.is_null(path));
        require(ar.dimensions(path) == 2);
        require(ar.extent(path) == std::vector<std::size_t>({2, 0}));
        ar.read(path, static_cast<int *>(nullptr), {2, 0});
    }
    for (auto path : {"/empty_string", "/group/@empty_string"}) {
        ar.write(path, static_cast<std::string const *>(nullptr), {2, 0});
        require(!ar.is_null(path));
        require(ar.extent(path) == std::vector<std::size_t>({2, 0}));
        ar.read(path, static_cast<std::string *>(nullptr), {2, 0});
    }
}

void invalid_selections(alps::hdf5::archive & ar) {
    std::array<int, 6> const values{1, 2, 3, 4, 5, 6};
    std::array<std::string, 6> const texts{"a", "b", "c", "d", "e", "f"};
    ar.write("/protected", values.data(), {2, 3});
    ar.write("/protected/@kept", 42);
    ar.write("/protected_group/child", 91);
    ar.write("/protected_group/@kept", 92);
    ar.write("/group/@protected", values.data(), {2, 3});
    ar.write("/group/@kept", 93);
    auto const preserved = [&] {
        std::array<int, 6> actual;
        ar.read("/protected", actual.data(), {2, 3});
        require(actual == values);
        ar.read("/group/@protected", actual.data(), {2, 3});
        require(actual == values);
        int scalar = 0;
        ar.read("/protected/@kept", scalar);
        require(scalar == 42);
        ar.read("/protected_group/child", scalar);
        require(scalar == 91 && ar.is_group("/protected_group"));
        ar.read("/protected_group/@kept", scalar);
        require(scalar == 92);
        ar.read("/group/@kept", scalar);
        require(scalar == 93);
    };
    struct selection {
        std::vector<std::size_t> size, chunk, offset;
    };
    auto const maximum = std::numeric_limits<std::size_t>::max();
    std::vector<selection> const invalid{
        {{2, 3}, {2}, {0, 0}},             // short chunk rank
        {{2, 3}, {2, 3, 1}, {0, 0}},       // long chunk rank
        {{2, 3}, {2, 3}, {0}},             // offset rank
        {{2, 3}, {1, 2}, {2, 2}},         // bounds
        {{2, 3}, {2, 3}, {1, 0}},         // full-sized chunk with nonzero offset
        {{2, 3}, {1, 1}, {maximum, 0}},   // offset + chunk overflow
        {{maximum, 2}, {1, 1}, {0, 0}},   // dataset element count overflow
        {{maximum, 2}, {maximum, 2}, {0, 0}} // transfer element count overflow
    };
    // These malformed requests must be rejected before allocation, conversion
    // or replacement. The buffers remain small even for overflowing extents.
    for (auto path : {"/protected", "/protected_group", "/group/@protected"}) {
        for (auto const & request : invalid) {
            rejects([&] { ar.write(path, values.data(), request.size, request.chunk, request.offset); });
            preserved();
            rejects([&] { ar.write(path, texts.data(), request.size, request.chunk, request.offset); });
            preserved();
        }
        rejects([&] { ar.write(path, static_cast<int const *>(nullptr), {2, 3}); });
        preserved();
        rejects([&] { ar.write(path, static_cast<std::string const *>(nullptr), {2, 3}); });
        preserved();
        rejects([&] {
            ar.write(path, values.data(), {maximum / sizeof(int) + 1, 1}, {1, 1}, {0, 0});
        });
        preserved();
    }
    rejects([&] { ar.write("/missing/@value", 1); });
    require(!ar.is_group("/missing"));
    rejects([&] { ar.write("/missing/@value", texts.data(), {2, 3}); });
    require(!ar.is_group("/missing"));
}
}


class Hdf5Write : public ::testing::Test {
protected:
    alps::testing::TemporaryDirectory temporary;
    std::string filename = (temporary.path() / "write.h5").string();
    alps::hdf5::archive ar{filename, "w"};
    void SetUp() override { ar.create_group("/group"); }
};

TEST_F(Hdf5Write, CompleteAndPartialTransfers) {
    array_writes<int>(ar, "numbers", {0, 1, 2, 3, 4, 5}, {7, 8});
    array_writes<double>(ar, "reals", {0.5, 1.5, 2.5, 3.5, 4.5, 5.5}, {7.5, 8.5});
    array_writes<bool>(ar, "flags", {false, true, false, true, false, true}, {true, true});
    array_writes<std::string>(ar, "strings", {"", "one", "two", "three", "four", "five"}, {"seven", "eight"});
}

TEST_F(Hdf5Write, ReplacementPreservesOnlyCompatibleMetadata) {
    ar.write("/numbers", 1);
    ar.write("/group/@sibling", 99);
    for (auto path : {"/group/@value", "/numbers/@value", "/@root", "@root"})
        attribute_replacements(ar, path);
    int sibling = 0;
    ar.read("/group/@sibling", sibling);
    EXPECT_EQ(sibling, 99);
    replacements(ar);
    ar.delete_attribute("/@root");
    EXPECT_FALSE(ar.is_attribute("@root"));
}

TEST_F(Hdf5Write, EmptyShapes) { empty_writes(ar); }

TEST_F(Hdf5Write, RejectedRequestsPreserveExistingData) { invalid_selections(ar); }

TEST_F(Hdf5Write, ClosedAndReadOnly) {
    std::array<int, 2> const values{1, 2};
    ar.write("/numbers", values.data(), {2});
    ar.close();
    EXPECT_THROW(ar.write("/closed", 1), alps::hdf5::archive_closed);
    EXPECT_THROW(ar.write("/closed", values.data(), {2}), alps::hdf5::archive_closed);
    alps::hdf5::archive input(filename, "r");
    EXPECT_THROW(input.write("/numbers", 1), alps::hdf5::archive_error);
    EXPECT_THROW(input.write("/numbers", values.data(), {2}), alps::hdf5::archive_error);
    std::array<int, 2> actual{};
    input.read("/numbers", actual.data(), {2});
    EXPECT_EQ(actual, values);
}
