// SPDX-License-Identifier: MIT
// Regression coverage for the logical path and segment-encoding contracts
// adapted from ALPS reconciliation 988651d41, src/alps/hdf5/src/archive.cpp.

#include <gtest/gtest.h>
#include <alps/hdf5/archive.hpp>
#include <alps/testing/temporary_directory.hpp>

#include <string>
#include <vector>

TEST(Hdf5Paths, ResolvesOnlyCompleteDotSegments) {
    alps::testing::TemporaryDirectory temporary;
    alps::hdf5::archive ar((temporary.path() / "paths.h5").string(), "a");
    ar.set_context("/outer/inner");
    EXPECT_EQ(ar.complete_path(""), "/outer/inner");
    EXPECT_EQ(ar.complete_path("."), "/outer/inner");
    EXPECT_EQ(ar.complete_path(".."), "/outer");
    EXPECT_EQ(ar.complete_path("../../../value"), "/value");
    EXPECT_EQ(ar.complete_path("child/../value"), "/outer/inner/value");
    EXPECT_EQ(ar.complete_path("..prefix"), "/outer/inner/..prefix");
    EXPECT_EQ(ar.complete_path(".../value"), "/outer/inner/.../value");
    EXPECT_EQ(ar.complete_path(".hidden"), "/outer/inner/.hidden");
    EXPECT_EQ(ar.complete_path("//absolute/./child/../value///"), "/absolute/value");
    EXPECT_EQ(ar.complete_path("/../../"), "/");
    EXPECT_EQ(ar.get_context(), "/outer/inner");
    ar.set_context("./child/../sibling//");
    EXPECT_EQ(ar.get_context(), "/outer/inner/sibling");
}

TEST(Hdf5Paths, DotPrefixedNamesRemainDistinctFromParentPaths) {
    alps::testing::TemporaryDirectory temporary;
    auto const filename = (temporary.path() / "dots.h5").string();
    {
        alps::hdf5::archive ar(filename, "a");
        ar.set_context("/outer/inner");
        ar.write("..value", 11);
        ar.write("../value", 22);
        ar.write("child/../value", 33);
        ar.write(".../value", 44);
        ar.close();
    }
    alps::hdf5::archive ar(filename, "r");
    int value = 0;
    ar.read("/outer/inner/..value", value);
    EXPECT_EQ(value, 11);
    ar.read("/outer/value", value);
    EXPECT_EQ(value, 22);
    ar.read("/outer/inner/value", value);
    EXPECT_EQ(value, 33);
    ar.read("/outer/inner/.../value", value);
    EXPECT_EQ(value, 44);
}

TEST(Hdf5Paths, PreservesAttributeNamesDuringParentNormalization) {
    alps::testing::TemporaryDirectory temporary;
    alps::hdf5::archive ar((temporary.path() / "attributes.h5").string(), "a");
    ar.create_group("/outer");
    ar.set_context("/outer/inner");
    EXPECT_EQ(ar.complete_path("../@part/../name/"), "/outer/@part/../name/");
    ar.write("../@part/../name/", 17);
    EXPECT_EQ(ar.list_attributes(".."), (std::vector<std::string>{"part/../name/"}));
    int value = 0;
    ar.read("/outer/./@part/../name/", value);
    EXPECT_EQ(value, 17);
    ar.write("/../@root", 23);
    ar.read("/@root", value);
    EXPECT_EQ(value, 23);
}

TEST(Hdf5Paths, EncodesExistingWireFormatAndDecodesOnce) {
    alps::testing::TemporaryDirectory temporary;
    alps::hdf5::archive ar((temporary.path() / "segments.h5").string(), "a");
    EXPECT_EQ(ar.encode_segment("A/B&C"), "A&#47;B&#38;C");
    EXPECT_EQ(ar.encode_segment("&#47;"), "&#38;#47;");
    EXPECT_EQ(ar.decode_segment("&#38;#47;"), "&#47;");
    EXPECT_EQ(ar.decode_segment("&#47;&#38;&#47;"), "/&/");
    std::vector<std::string> names{
        "A/B&C", "&#47;", "&#38;", "&&//", "a&b", "&#", "&#not-a-number;",
        "&#65;", "&#47", "&#038;", "&#999999999999999999999999;", "plain"
    };
    int expected = 0;
    for (auto const & name : names) {
        auto const encoded = ar.encode_segment(name);
        EXPECT_EQ(ar.decode_segment(encoded), name);
        ar.write("/names/" + encoded, ++expected);
    }
    expected = 0;
    for (auto const & name : names) {
        int value = 0;
        ar.read("/names/" + ar.encode_segment(name), value);
        EXPECT_EQ(value, ++expected);
    }
    EXPECT_EQ(ar.list_children("/names").size(), names.size());
}

TEST(Hdf5Paths, LeavesUnknownAndMalformedEscapesLiteral) {
    alps::testing::TemporaryDirectory temporary;
    alps::hdf5::archive ar((temporary.path() / "literal.h5").string(), "a");
    for (std::string const text : {"&", "a&b", "&name;", "&#", "&#not-a-number;",
                                   "&#65;", "&#47", "&#038;", "&#-1;",
                                   "&#999999999999999999999999;"})
        EXPECT_EQ(ar.decode_segment(text), text);
    EXPECT_EQ(ar.decode_segment("&literal/&#47;&#unknown;&#38;"), "&literal//&#unknown;&");
}
