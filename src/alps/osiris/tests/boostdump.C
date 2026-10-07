/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2003 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <gtest/gtest.h>
#include <alps/osiris.h>
#include <alps/testing/temporary_directory.hpp>
#include <alps/osiris/archivedump.h>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <fstream>

TEST(OsirisSerialization, TextArchivePreservesScalarValues) {
  alps::testing::TemporaryDirectory directory;
  const auto path = directory.path() / "record.dump";
  {
    std::ofstream stream(path);
    ASSERT_TRUE(stream.good());
    boost::archive::text_oarchive archive(stream);
    alps::archive_odump<boost::archive::text_oarchive> output(archive);
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
  std::ifstream stream(path);
  ASSERT_TRUE(stream.good());
  boost::archive::text_iarchive archive(stream);
  alps::archive_idump<boost::archive::text_iarchive> id(archive);

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
}
