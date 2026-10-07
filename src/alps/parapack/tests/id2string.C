/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2008 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <alps/parapack/util.h>
#include <gtest/gtest.h>
#include <string>
#include <utility>

namespace {
class ProcessIdentifier : public ::testing::TestWithParam<std::pair<int, std::string>> {};
TEST_P(ProcessIdentifier, KeepsHistoricalSortableSpelling) {
  EXPECT_EQ(alps::id2string(GetParam().first), GetParam().second);
}
INSTANTIATE_TEST_SUITE_P(
    HistoricalValues, ProcessIdentifier,
    ::testing::Values(std::make_pair(0, "0"), std::make_pair(9, "9"), std::make_pair(10, "_10"),
                      std::make_pair(15, "_15"), std::make_pair(100, "__100"),
                      std::make_pair(200, "__200"), std::make_pair(1000, "___1000"),
                      std::make_pair(1001, "___1001"), std::make_pair(100000, "_____100000")));
} // namespace
