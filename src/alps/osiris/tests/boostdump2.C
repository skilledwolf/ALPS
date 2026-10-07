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

TEST(OsirisSerialization, ArchiveAdapterReadsHistoricalXdrRecord) {
  // Immutable legacy wire fixture: this check is independent of today's writer.
  alps::IXDRFileDump id(boost::filesystem::path(ALPS_TEST_SOURCE_DIR "/xdrdump2.dump"));
  alps::idump_archive archive(id);

  EXPECT_EQ(alps::get<bool>(archive), false);
  EXPECT_EQ(alps::get<int8_t>(archive), 63);
  EXPECT_EQ(alps::get<uint8_t>(archive), 201);
  EXPECT_EQ(alps::get<int16_t>(archive), -699);
  EXPECT_EQ(alps::get<uint16_t>(archive), 43299);
  EXPECT_EQ(alps::get<int32_t>(archive), 847229);
  EXPECT_EQ(alps::get<uint32_t>(archive), 4294967295u);
  EXPECT_EQ(alps::get<int64_t>(archive), -1152921504606846976ll);
  EXPECT_EQ(alps::get<uint64_t>(archive), 18446744073709551614ull);
  // The wire and text archives preserve the full double, not six printed digits.
  EXPECT_EQ(alps::get<double>(archive), 3.14159265358979323846);
  std::string text;
  archive >> text;
  EXPECT_EQ(text, "test string");
}
