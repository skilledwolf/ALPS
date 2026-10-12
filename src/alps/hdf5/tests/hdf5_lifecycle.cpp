// SPDX-License-Identifier: MIT
// Adapted from the ownership/finalization contracts of ALPS reconciliation
// 81e7ad00d. These tests retain the current SDK's per-view close semantics.
#include <alps/hdf5/archive.hpp>
#include <alps/testing/temporary_directory.hpp>
#include <gtest/gtest.h>
#include <hdf5.h>
#include <boost/filesystem.hpp>
#include <memory>

class ArchiveLifetime : public ::testing::Test {
protected:
    alps::testing::TemporaryDirectory temporary;
    std::string filename = (temporary.path() / "lifetime.h5").string();
};

TEST_F(ArchiveLifetime, CopiesKeepFileAliveWithIndependentContexts) {
    alps::hdf5::archive owner(filename, "a");
    owner.write("/first/value", 17);
    owner.write("/second/value", 23);
    owner.set_context("/first");
    alps::hdf5::archive copy(owner);
    copy.set_context("/second");
    EXPECT_EQ(owner.get_context(), "/first");
    owner.close();
    EXPECT_FALSE(owner.is_open());
    ASSERT_TRUE(copy.is_open());
    int value = 0;
    copy.read("value", value);
    EXPECT_EQ(value, 23);
    copy.write("new", 31);
    copy.close();
    alps::hdf5::archive closed_copy(owner);
    EXPECT_FALSE(closed_copy.is_open());
    EXPECT_THROW(closed_copy.read("/first/value", value), alps::hdf5::archive_closed);
    alps::hdf5::archive reopened(filename, "r");
    reopened.read("/second/new", value);
    EXPECT_EQ(value, 31);
}

TEST_F(ArchiveLifetime, ProxyRetainsItsView) {
    alps::hdf5::archive owner(filename, "a");
    auto proxy = owner["/value"];
    owner.close();
    proxy << 29;
    int value = 0;
    proxy >> value;
    EXPECT_EQ(value, 29);
}

TEST_F(ArchiveLifetime, ReopeningForWriteUpgradesTheSharedContext) {
    { alps::hdf5::archive seed(filename, "a"); seed.write("/value", 1); }
    alps::hdf5::archive first(filename, "r");
    {
        alps::hdf5::archive writer(filename, "a");
        writer.write("/value", 7);
        writer.close();
    }
    int value = 0;
    first.read("/value", value);
    EXPECT_EQ(value, 7);
    first.close();
}

TEST_F(ArchiveLifetime, ExplicitCloseReportsOutstandingHdf5Objects) {
    alps::hdf5::archive ar(filename, "a");
    ar.write("/value", 7);
    // Borrow the sole archive file id without adding another file reference.
    hid_t file = -1;
    ASSERT_EQ(H5Fget_obj_ids(H5F_OBJ_ALL, H5F_OBJ_FILE, 1, &file), 1);
    hid_t dataset = H5Dopen2(file, "/value", H5P_DEFAULT);
    ASSERT_GE(dataset, 0);
    // H5F_CLOSE_SEMI refuses the last file close while a dataset remains open.
    EXPECT_THROW(ar.close(), alps::hdf5::archive_error);
    EXPECT_TRUE(ar.is_open());
    ASSERT_GE(H5Dclose(dataset), 0);
    EXPECT_NO_THROW(ar.close());
    EXPECT_FALSE(ar.is_open());
}

TEST_F(ArchiveLifetime, FailedReplacementLeavesDestinationIntact) {
    alps::hdf5::archive ar(filename, "w");
    ar.write("/value", 7);
    // Force the publication rename to fail without needing a full disk or
    // elevated permissions. A nonempty directory is a protected destination.
    boost::filesystem::create_directory(filename);
    auto const marker = boost::filesystem::path(filename) / "keep";
    boost::filesystem::create_directory(marker);
    EXPECT_THROW(ar.close(), boost::filesystem::filesystem_error);
    EXPECT_TRUE(boost::filesystem::is_directory(marker));
    EXPECT_FALSE(ar.is_open());
}

TEST_F(ArchiveLifetime, RelativeAndAbsoluteNamesShareOneOwner) {
    { alps::hdf5::archive seed(filename, "a"); seed.write("/value", 1); }
    auto const relative = boost::filesystem::relative(filename, boost::filesystem::current_path());
    alps::hdf5::archive reader(relative, "r");
    EXPECT_EQ(reader.get_filename(), relative.string());
    alps::hdf5::archive writer(filename, "a");
    writer.write("/value", 17);
    writer.close();
    int value = 0;
    reader.read("/value", value);
    EXPECT_EQ(value, 17);
}
