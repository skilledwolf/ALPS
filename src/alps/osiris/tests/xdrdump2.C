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

TEST(OsirisSerialization, ReadsHistoricalXdrRecord) {
  // Immutable legacy wire fixture: this check is independent of today's writer.
  alps::IXDRFileDump id(boost::filesystem::path(ALPS_TEST_SOURCE_DIR "/xdrdump2.dump"));

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
