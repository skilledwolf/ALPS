// SPDX-License-Identifier: MIT
// Adapted from ALPS reconciliation 81e7ad00d, mc/tests/checkpoint.cpp.
#include <alps/ngs.hpp>
#include <alps/testing/temporary_directory.hpp>
#include <gtest/gtest.h>

namespace {
alps::mcobservables observations() {
    alps::mcobservables values;
    values << alps::accumulator::RealObservable("kept") << alps::accumulator::RealObservable("gone");
    values["kept"] << 1.;
    values["gone"] << 2.;
    return values;
}
alps::mcresults results() {
    auto measured = observations();
    alps::mcresults values;
    values.insert("kept", alps::mcresult(measured["kept"]));
    values.insert("gone", alps::mcresult(measured["gone"]));
    return values;
}

template<class Collection> void replacement(Collection values) {
    alps::testing::TemporaryDirectory temporary;
    alps::hdf5::archive ar((temporary.path() / "collections.h5").string(), "a");
    ar["/unrelated"] << 8;
    ar["/collection"] << values;
    ar["/collection/@tag"] << 9;
    ar["/collection/orphan"] << 10;
    values.erase("gone");
    ar["/collection"] << values;
    Collection restored;
    ar["/collection"] >> restored;
    ASSERT_EQ(restored.size(), 1u);
    EXPECT_TRUE(restored.has("kept"));
    EXPECT_FALSE(restored.has("gone"));
    EXPECT_FALSE(ar.is_data("/collection/orphan"));
    int value = 0;
    ar["/unrelated"] >> value;
    EXPECT_EQ(value, 8);
    ar["/collection/@tag"] >> value;
    EXPECT_EQ(value, 9);
    values.clear();
    ar["/collection"] << values;
    Collection empty;
    ar["/collection"] >> empty;
    EXPECT_TRUE(empty.empty());
    EXPECT_TRUE(ar.is_group("/collection"));
    EXPECT_TRUE(ar.list_children("/collection").empty());
    ar["/collection/@tag"] >> value;
    EXPECT_EQ(value, 9);
    // A caller choosing root as the collection context owns all root children.
    ar["/"] << restored;
    ar["/@tag"] << 11;
    ar["/"] << values;
    EXPECT_TRUE(ar.list_children("/").empty());
    ar["/@tag"] >> value;
    EXPECT_EQ(value, 11);
}
}

TEST(CollectionCheckpoint, ObservablesReplaceOwnedChildren) { replacement(observations()); }
TEST(CollectionCheckpoint, ResultsReplaceOwnedChildren) { replacement(results()); }

TEST(CollectionCheckpoint, SavingEmptyResultsUpdatesOnlyTheRequestedGroup) {
    alps::testing::TemporaryDirectory temporary;
    auto const filename = (temporary.path() / "results.h5").string();
    { alps::hdf5::archive ar(filename, "a"); ar["/unrelated"] << 8; }
    alps::params parameters;
    parameters["revision"] = 1;
    auto const values = results();
    alps::save_results(values, parameters, filename, "/first/results");
    alps::save_results(values, parameters, filename, "/second/results");
    parameters["revision"] = 2;
    alps::save_results(alps::mcresults{}, parameters, filename, "/first/results");
    alps::hdf5::archive ar(filename);
    alps::mcresults empty, retained;
    ar["/first/results"] >> empty;
    ar["/second/results"] >> retained;
    EXPECT_TRUE(empty.empty());
    EXPECT_TRUE(ar.is_group("/first/results"));
    EXPECT_TRUE(ar.list_children("/first/results").empty());
    ASSERT_EQ(retained.size(), 2u);
    EXPECT_EQ(retained["kept"].count(), 1u);
    alps::params stored(ar, "/parameters");
    EXPECT_EQ(stored["revision"].cast<int>(), 2);
    int unrelated = 0;
    ar["/unrelated"] >> unrelated;
    EXPECT_EQ(unrelated, 8);
}

TEST(CollectionCheckpoint, ResultPathsCannotEraseParameterMetadata) {
    alps::testing::TemporaryDirectory temporary;
    auto const filename = (temporary.path() / "protected.h5").string();
    alps::params old;
    old["kept"] = 17;
    { alps::hdf5::archive ar(filename, "a"); ar["/parameters"] << old; }
    alps::params changed;
    changed["changed"] = 23;
    for (auto path : {"/", "/parameters", "/parameters/child"}) {
        EXPECT_THROW(alps::save_results(results(), changed, filename, path), alps::hdf5::invalid_path);
        alps::hdf5::archive ar(filename);
        alps::params preserved(ar, "/parameters");
        ASSERT_EQ(preserved.size(), 1u);
        EXPECT_EQ(preserved["kept"].cast<int>(), 17);
    }
}
