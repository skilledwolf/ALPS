/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2005 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <gtest/gtest.h>
#include <alps/osiris.h>
#include <alps/testing/temporary_directory.hpp>
#include <fstream>
#include <iterator>

TEST(OsirisSerialization, NativeDumpReadsArchiveAdapter) {
  alps::testing::TemporaryDirectory directory;
  const auto path = directory.path() / "record.dump";
  {
    alps::OXDRFileDump dump(boost::filesystem::path(path.string()));
    alps::odump_archive output(dump);
    output << bool(false);
    output << int8_t(63);
    output << uint8_t(201);
    output << int16_t(-699);
    output << uint16_t(43299);
    output << int32_t(847229);
    output << uint32_t(4294967295u);
    output << int64_t(-1152921504606846976ll);
    output << uint64_t(18446744073709551614ull);
    output << double(3.14159265358979323846);
    output << std::string("test string");
  }
  alps::IXDRFileDump id(boost::filesystem::path(path.string()));

  EXPECT_EQ(id.get<bool>(), false);
  EXPECT_EQ(id.get<int8_t>(), 63);
  EXPECT_EQ(id.get<uint8_t>(), 201);
  EXPECT_EQ(id.get<int16_t>(), -699);
  EXPECT_EQ(id.get<uint16_t>(), 43299);
  EXPECT_EQ(static_cast<int32_t>(id), 847229);
  EXPECT_EQ(static_cast<uint32_t>(id), 4294967295u);
  EXPECT_EQ(static_cast<int64_t>(id), -1152921504606846976ll);
  EXPECT_EQ(static_cast<uint64_t>(id), 18446744073709551614ull);
  // The wire and text archives preserve the full double, not six printed digits.
  EXPECT_EQ(static_cast<double>(id), 3.14159265358979323846);
  std::string text;
  id >> text;
  EXPECT_EQ(text, "test string");

  // A matching writer/reader bug must not silently change the checkpoint format.
  std::ifstream expected(ALPS_TEST_SOURCE_DIR "/xdrdump2.dump", std::ios::binary);
  std::ifstream actual(path, std::ios::binary);
  ASSERT_TRUE(expected.good());
  ASSERT_TRUE(actual.good());
  const std::string expected_bytes{std::istreambuf_iterator<char>(expected), {}};
  const std::string actual_bytes{std::istreambuf_iterator<char>(actual), {}};
  EXPECT_EQ(actual_bytes, expected_bytes);
}
