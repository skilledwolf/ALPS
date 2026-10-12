// SPDX-License-Identifier: MIT
// Regression coverage for the sequence-adapter consolidation adapted from
// ALPS reconciliation 2cc35a558, including the retained ALPS storage contracts.

#include <gtest/gtest.h>
#include <alps/hdf5.hpp>
#include <alps/hdf5/array.hpp>
#include <alps/hdf5/stdarray.hpp>
#include <alps/testing/temporary_directory.hpp>

#include <array>
#include <complex>
#include <limits>
#include <string>
#include <valarray>
#include <vector>

namespace {
template<class C> void expect_values(C const& actual, std::vector<typename C::value_type> const& expected) {
    ASSERT_EQ(actual.size(), expected.size());
    for (std::size_t i = 0; i < actual.size(); ++i) EXPECT_EQ(actual[i], expected[i]);
}
template<class C> void round_trip(alps::hdf5::archive& ar, std::string const& path,
                                  C const& input, C output) {
    ar[path] = input;
    ar[path] >> output;
    ASSERT_EQ(output.size(), input.size());
    for (std::size_t i = 0; i < input.size(); ++i) EXPECT_EQ(output[i], input[i]);
}
}

TEST(Hdf5Sequence, RoundTripsEachContainerAndRetainsComplexMarkers) {
    alps::testing::TemporaryDirectory temporary;
    alps::hdf5::archive ar((temporary.path() / "roundtrip.h5").string(), "a");
    round_trip(ar, "/vector", std::vector<int>{2, 3, 5}, std::vector<int>{});
    round_trip(ar, "/stdarray", std::array<int, 3>{{2, 3, 5}}, std::array<int, 3>{});
    round_trip(ar, "/boostarray", boost::array<int, 3>{{2, 3, 5}}, boost::array<int, 3>{});
    round_trip(ar, "/valarray", std::valarray<int>{2, 3, 5}, std::valarray<int>{});
    boost::numeric::ublas::vector<int> ublas(3);
    ublas[0] = 2; ublas[1] = 3; ublas[2] = 5;
    round_trip(ar, "/ublas", ublas, boost::numeric::ublas::vector<int>{});
    using complex = std::complex<double>;
    const complex first(1, -2), second(3, 4);
    round_trip(ar, "/complex_vector", std::vector<complex>{first, second}, std::vector<complex>{});
    round_trip(ar, "/complex_stdarray", std::array<complex, 2>{{first, second}}, std::array<complex, 2>{});
    round_trip(ar, "/complex_boostarray", boost::array<complex, 2>{{first, second}}, boost::array<complex, 2>{});
    round_trip(ar, "/complex_valarray", std::valarray<complex>{first, second}, std::valarray<complex>{});
    boost::numeric::ublas::vector<complex> complex_ublas(2);
    complex_ublas[0] = first; complex_ublas[1] = second;
    round_trip(ar, "/complex_ublas", complex_ublas, boost::numeric::ublas::vector<complex>{});
    for (auto const name : {"vector", "stdarray", "boostarray", "valarray", "ublas"})
        EXPECT_TRUE(ar.is_complex(std::string("/complex_") + name));
    std::vector<double> real{17};
    EXPECT_THROW(ar["/complex_vector"] >> real, alps::hdf5::archive_error);
    EXPECT_EQ(real, (std::vector<double>{17}));
    boost::numeric::ublas::vector<double> real_ublas(1);
    real_ublas[0] = 17;
    EXPECT_THROW(ar["/complex_ublas"] >> real_ublas, alps::hdf5::archive_error);
    expect_values(real_ublas, {17});
    std::vector<complex> complex_output{first};
    EXPECT_THROW(ar["/vector"] >> complex_output, alps::hdf5::archive_error);
    EXPECT_EQ(complex_output, (std::vector<complex>{first}));
}

TEST(Hdf5Sequence, KeepsLegacyTopLevelEmptyNullValues) {
    alps::testing::TemporaryDirectory temporary;
    alps::hdf5::archive ar((temporary.path() / "empty.h5").string(), "a");
    round_trip(ar, "/vector", std::vector<int>{}, std::vector<int>{17});
    round_trip(ar, "/stdarray", std::array<int, 0>{}, std::array<int, 0>{});
    round_trip(ar, "/boostarray", boost::array<int, 0>{}, boost::array<int, 0>{});
    round_trip(ar, "/valarray", std::valarray<int>{}, std::valarray<int>{17});
    round_trip(ar, "/ublas", boost::numeric::ublas::vector<int>{}, boost::numeric::ublas::vector<int>(1));
    round_trip(ar, "/bool", std::vector<bool>{}, std::vector<bool>{true});
    std::vector<std::vector<int>> nested, nested_output{{17}};
    ar["/nested"] = nested;
    ar["/nested"] >> nested_output;
    EXPECT_TRUE(nested_output.empty());
    using complex = std::complex<double>;
    round_trip(ar, "/complex", std::vector<complex>{}, std::vector<complex>{{1, 2}});
    EXPECT_TRUE(ar.is_complex("/complex"));
    for (auto const name : {"vector", "stdarray", "boostarray", "valarray", "ublas", "bool", "nested", "complex"})
        EXPECT_TRUE(ar.is_null(std::string("/") + name));
    EXPECT_TRUE(ar.is_datatype<int>("/nested"));
    EXPECT_EQ(alps::hdf5::get_pointer(std::vector<int>{}), nullptr);
    const std::array<int, 0> empty_array{};
    EXPECT_EQ(alps::hdf5::get_pointer(empty_array), nullptr);
    EXPECT_TRUE(alps::hdf5::is_vectorizable(std::array<std::vector<int>, 0>{}));
}

TEST(Hdf5Sequence, PreservesRankedEmptyInnerSequences) {
    alps::testing::TemporaryDirectory temporary;
    alps::hdf5::archive ar((temporary.path() / "ranked-empty.h5").string(), "a");
    std::vector<std::vector<int>> input(3), output;
    ar["/vector"] = input;
    EXPECT_FALSE(ar.is_null("/vector"));
    EXPECT_EQ(ar.extent("/vector"), (std::vector<std::size_t>{3, 0}));
    ar["/vector"] >> output;
    EXPECT_EQ(output, input);
    const std::array<std::array<int, 0>, 3> fixed{};
    std::array<std::array<int, 0>, 3> fixed_output{};
    ar["/fixed"] = fixed;
    EXPECT_EQ(ar.extent("/fixed"), (std::vector<std::size_t>{3, 0}));
    ar["/fixed"] >> fixed_output;
    EXPECT_EQ(fixed_output, fixed);
    std::vector<std::vector<bool>> booleans(3), restored_booleans;
    ar["/booleans"] = booleans;
    ar["/booleans"] >> restored_booleans;
    EXPECT_EQ(ar.extent("/booleans"), (std::vector<std::size_t>{3, 0}));
    EXPECT_EQ(restored_booleans, booleans);
    std::vector<std::vector<std::complex<double>>> complex_rows(3), restored_complex;
    ar["/complex"] = complex_rows;
    ar["/complex"] >> restored_complex;
    EXPECT_EQ(ar.extent("/complex"), (std::vector<std::size_t>{3, 0, 2}));
    EXPECT_EQ(restored_complex, complex_rows);
}

TEST(Hdf5Sequence, BulkBooleanTransfersSupportAttributesAndNestedRows) {
    alps::testing::TemporaryDirectory temporary;
    alps::hdf5::archive ar((temporary.path() / "booleans.h5").string(), "a");
    std::vector<bool> input(513);
    for (std::size_t i = 0; i < input.size(); ++i) input[i] = i % 3 == 0;
    round_trip(ar, "/value", input, std::vector<bool>{});
    round_trip(ar, "/value/@flags", input, std::vector<bool>{});
    EXPECT_EQ(ar.extent("/value/@flags"), (std::vector<std::size_t>{513}));
    EXPECT_FALSE(ar.is_complex("/value/@flags"));
    std::vector<std::vector<bool>> rows{{true, false, true}, {false, true, false}}, output;
    ar["/rows"] = rows;
    ar["/rows"] >> output;
    EXPECT_EQ(output, rows);
}

TEST(Hdf5Sequence, RejectsWrongRankAndFixedSizesBeforeChangingTheContainer) {
    alps::testing::TemporaryDirectory temporary;
    alps::hdf5::archive ar((temporary.path() / "dimensions.h5").string(), "a");
    ar["/scalar"] = 3;
    ar["/three"] = std::vector<int>{1, 2, 3};
    ar["/matrix"] = std::vector<std::vector<int>>{{1, 2, 3}, {4, 5, 6}};
    std::vector<int> vector{17};
    EXPECT_THROW(ar["/scalar"] >> vector, alps::hdf5::archive_error);
    EXPECT_THROW(ar["/matrix"] >> vector, alps::hdf5::archive_error);
    EXPECT_EQ(vector, (std::vector<int>{17}));
    std::array<int, 2> fixed{{17, 23}};
    EXPECT_THROW(ar["/three"] >> fixed, alps::hdf5::archive_error);
    EXPECT_EQ(fixed, (std::array<int, 2>{{17, 23}}));
    boost::array<int, 2> boost_fixed{{17, 23}};
    EXPECT_THROW(ar["/three"] >> boost_fixed, alps::hdf5::archive_error);
    expect_values(boost_fixed, {17, 23});
    std::vector<std::array<int, 2>> inner_fixed{{{17, 23}}};
    EXPECT_THROW(ar["/matrix"] >> inner_fixed, alps::hdf5::archive_error);
    ASSERT_EQ(inner_fixed.size(), 1);
    EXPECT_EQ(inner_fixed.front(), fixed);
    std::vector<bool> booleans{true};
    EXPECT_THROW(ar["/matrix"] >> booleans, alps::hdf5::archive_error);
    EXPECT_EQ(booleans, (std::vector<bool>{true}));
    std::valarray<int> valarray{17};
    EXPECT_THROW(ar["/scalar"] >> valarray, alps::hdf5::archive_error);
    expect_values(valarray, {17});
    boost::numeric::ublas::vector<int> ublas(1);
    ublas[0] = 17;
    EXPECT_THROW(ar["/scalar"] >> ublas, alps::hdf5::archive_error);
    expect_values(ublas, {17});
    ar["/empty"] = std::vector<int>{};
    EXPECT_THROW(ar["/empty"] >> fixed, alps::hdf5::archive_error);
    EXPECT_EQ(fixed, (std::array<int, 2>{{17, 23}}));
}

TEST(Hdf5Sequence, RejectsMalformedSelectionsBeforeIndexingOrReplacing) {
    alps::testing::TemporaryDirectory temporary;
    alps::hdf5::archive ar((temporary.path() / "selections.h5").string(), "a");
    ar["/matrix"] = std::vector<std::vector<int>>{{1, 2, 3}, {4, 5, 6}};
    std::vector<int> output{17};
    EXPECT_THROW(alps::hdf5::load(ar, "/matrix", output, {1, 1}, {0, 0}), alps::hdf5::archive_error);
    EXPECT_THROW(alps::hdf5::load(ar, "/matrix", output, {1}, {}), alps::hdf5::archive_error);
    EXPECT_THROW(alps::hdf5::load(ar, "/matrix", output, {}, {0, 0, 0}), alps::hdf5::archive_error);
    EXPECT_THROW(alps::hdf5::load(ar, "/matrix", output, {2}, {0}), alps::hdf5::archive_error);
    EXPECT_THROW(alps::hdf5::load(ar, "/matrix", output, {1}, {(std::numeric_limits<std::size_t>::max)()}), alps::hdf5::archive_error);
    EXPECT_EQ(output, (std::vector<int>{17}));
    alps::hdf5::load(ar, "/matrix", output, {1}, {1});
    EXPECT_EQ(output, (std::vector<int>{4, 5, 6}));
    ar["/old/sentinel"] = 23;
    EXPECT_THROW(alps::hdf5::save(ar, "/old", output, {2}, {}, {}), alps::hdf5::archive_error);
    EXPECT_THROW(alps::hdf5::save(ar, "/old", output, {2}, {2}, {0}), alps::hdf5::archive_error);
    EXPECT_THROW(alps::hdf5::save(ar, "/old", output, {2}, {1}, {2}), alps::hdf5::archive_error);
    int sentinel = 0;
    ar["/old/sentinel"] >> sentinel;
    EXPECT_EQ(sentinel, 23);
}

TEST(Hdf5Sequence, RejectsNoncanonicalGroupIndicesBeforeResizing) {
    alps::testing::TemporaryDirectory temporary;
    alps::hdf5::archive ar((temporary.path() / "indices.h5").string(), "a");
    for (std::string const invalid : {"2", "01", "-1", "+1", "name", "999999999999999999999999"}) {
        SCOPED_TRACE(invalid);
        if (ar.is_group("/values")) ar.delete_group("/values");
        ar["/values/0"] = 2;
        ar["/values/" + invalid] = 3;
        std::vector<int> output{17};
        EXPECT_THROW(ar["/values"] >> output, alps::hdf5::invalid_path);
        EXPECT_EQ(output, (std::vector<int>{17}));
        std::array<int, 2> fixed{{17, 23}};
        EXPECT_THROW(ar["/values"] >> fixed, alps::hdf5::invalid_path);
        EXPECT_EQ(fixed, (std::array<int, 2>{{17, 23}}));
    }
}

TEST(Hdf5Sequence, RestoresRaggedGroupsInNumericOrderAndRemovesOldMembers) {
    alps::testing::TemporaryDirectory temporary;
    alps::hdf5::archive ar((temporary.path() / "ragged.h5").string(), "a");
    std::vector<std::vector<int>> input(12), output;
    for (std::size_t i = 0; i < input.size(); ++i) input[i].assign(i + 1, static_cast<int>(i));
    ar["/value"] = input;
    ASSERT_TRUE(ar.is_group("/value"));
    ar["/value"] >> output;
    EXPECT_EQ(output, input);
    input.resize(2);
    ar["/value"] = input;
    EXPECT_EQ(ar.list_children("/value"), (std::vector<std::string>{"0", "1"}));
    ar["/value"] >> output;
    EXPECT_EQ(output, input);
}

TEST(Hdf5Sequence, ChecksNestedValarrayShapesWithoutComparingPastAnExtent) {
    std::valarray<std::vector<std::vector<int>>> values(2);
    values[0] = {{1}};
    values[1] = {};
    EXPECT_FALSE(alps::hdf5::is_vectorizable(values));
    EXPECT_THROW(alps::hdf5::get_extent(values), alps::hdf5::archive_error);
    std::valarray<int> empty;
    EXPECT_EQ(alps::hdf5::get_pointer(empty), nullptr);
    EXPECT_THROW(alps::hdf5::set_extent(empty, {}), alps::hdf5::archive_error);
}
