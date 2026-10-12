// SPDX-License-Identifier: MIT
// Adapted from ALPS reconciliation 81e7ad00d, mc/tests/checkpoint.cpp.
#include <alps/ngs.hpp>
#include <alps/testing/temporary_directory.hpp>
#include <gtest/gtest.h>
#include <hdf5.h>

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

namespace {
alps::params parameters() {
    alps::params values;
    values["kept"] = 1;
    values["gone"] = 2;
    return values;
}

template<class Collection> void rejects_native_children(Collection const &values) {
    alps::testing::TemporaryDirectory temporary;
    for (std::string const kind : {"dangling", "type", "zz@attribute", ".."}) {
        SCOPED_TRACE(kind);
        const auto filename = (temporary.path() / (kind + ".h5")).string();
        const auto child = "/collection/" + (kind == "type" || kind == "dangling"
                                              ? std::string("zz_unsupported") : kind);
        {
            alps::hdf5::archive ar(filename, "a");
            ar["/collection"] << values;
            ar["/collection/a_preserved"] << 17;
            ar["/collection/@owner"] << 23;
            ar["/sibling"] << 31;
        }
        auto file = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
        ASSERT_GE(file, 0);
        if (kind == "dangling") {
            ASSERT_GE(H5Lcreate_soft("/missing", file, child.c_str(), H5P_DEFAULT, H5P_DEFAULT), 0);
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
            ASSERT_GE(H5Dclose(data), 0);
            ASSERT_GE(H5Sclose(space), 0);
        }
        ASSERT_GE(H5Fclose(file), 0);
        {
            alps::hdf5::archive ar(filename, "a");
            EXPECT_THROW(ar["/collection"] << Collection{}, alps::hdf5::invalid_path);
            int value = 0;
            ar["/collection/a_preserved"] >> value; EXPECT_EQ(value, 17);
            ar["/collection/@owner"] >> value; EXPECT_EQ(value, 23);
            ar["/sibling"] >> value; EXPECT_EQ(value, 31);
            EXPECT_TRUE(ar.is_group("/collection/kept") || ar.is_data("/collection/kept"));
            EXPECT_TRUE(ar.is_group("/collection/gone") || ar.is_data("/collection/gone"));
            auto loaded = values;
            EXPECT_THROW(ar["/collection"] >> loaded, alps::hdf5::invalid_path);
            EXPECT_EQ(loaded.size(), values.size());
        }
        file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        ASSERT_GE(file, 0);
        EXPECT_GT(H5Lexists(file, child.c_str(), H5P_DEFAULT), 0);
        ASSERT_GE(H5Fclose(file), 0);
    }
}

template<class Collection> void rejects_dataset_ancestor() {
    alps::testing::TemporaryDirectory temporary;
    alps::hdf5::archive ar((temporary.path() / "ancestor.h5").string(), "a");
    ar["/ancestor"] << 17;
    ar["/ancestor/@owner"] << 23;
    EXPECT_THROW(ar["/ancestor/nested/collection"] << Collection{}, alps::hdf5::invalid_path);
    int value = 0;
    ar["/ancestor"] >> value; EXPECT_EQ(value, 17);
    ar["/ancestor/@owner"] >> value; EXPECT_EQ(value, 23);
    EXPECT_FALSE(ar.is_group("/ancestor"));
}

template<class Collection> void rejects_incoming_names(Collection const &invalid) {
    alps::testing::TemporaryDirectory temporary;
    alps::hdf5::archive ar((temporary.path() / "names.h5").string(), "a");
    ar["/collection/keep"] << 17;
    EXPECT_THROW(ar["/collection"] << invalid, alps::hdf5::invalid_path);
    int value = 0;
    ar["/collection/keep"] >> value;
    EXPECT_EQ(value, 17);
}
}

TEST(CollectionCheckpoint, NativeChildrenPreserveObservablesOnFailure) {
    rejects_native_children(observations());
    rejects_dataset_ancestor<alps::mcobservables>();
}

TEST(CollectionCheckpoint, NativeChildrenPreserveResultsOnFailure) {
    rejects_native_children(results());
    rejects_dataset_ancestor<alps::mcresults>();
}

TEST(CollectionCheckpoint, NativeChildrenPreserveParametersOnFailure) {
    rejects_native_children(parameters());
    rejects_dataset_ancestor<alps::params>();
}

TEST(CollectionCheckpoint, InvalidIncomingNamesPreservePreviousChildren) {
    for (const auto &name : std::vector<std::string>{"", ".", "..", "a@b", std::string("a\0b", 3)}) {
        SCOPED_TRACE(name);
        alps::mcobservables invalid;
        invalid << alps::accumulator::RealObservable(name);
        invalid[name] << 1.;
        rejects_incoming_names(invalid);
        alps::mcresults invalid_results;
        invalid_results.insert(name, alps::mcresult(invalid[name]));
        rejects_incoming_names(invalid_results);
    }
}
