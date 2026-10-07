#include <alps/testing/stream_fixture.hpp>
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

/* $Id$ */

#include <alps/osiris.h>
#include <iostream>
#include <cstdlib>

TEST(OsirisSerialization, Boostdump2) {
    alps::testing::StreamFixture transcript(ALPS_TEST_BINARY_DIR "/boostdump2.input");
    { // Flush serialization objects before checking the captured stream.

#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  std::string file;
  std::cin >> file;

  alps::IXDRFileDump id=alps::IXDRFileDump(boost::filesystem::path(file));
  alps::idump_archive ar(id);
  std::cout << alps::get<bool>(ar) << ' ';
  std::cout << static_cast<int32_t>(alps::get<int8_t>(ar)) << ' ';
  std::cout << static_cast<int32_t>(alps::get<uint8_t>(ar)) << ' ';
  std::cout << alps::get<int16_t>(ar) << ' ';
  std::cout << alps::get<uint16_t>(ar) << ' ';
  std::cout << alps::get<int32_t>(ar) << ' ';
  std::cout << alps::get<uint32_t>(ar) << ' ';
  int64_t i8 = alps::get<int64_t>(ar);
  uint64_t i9(alps::get<uint64_t>(ar));
  std::cout << i8 << ' '  << i9  << ' ';
  double i10 = alps::get<double>(ar);
  std::cout << i10 << ' ';
  std::string str;
  ar >> str;
  std::cout << str << std::endl;
  
#ifndef BOOST_NO_EXCEPTIONS
}
catch (std::exception& exp) {
  std::cerr << exp.what() << std::endl;
  FAIL() << "Unexpected exception in serialization contract";
}
#endif
    }
  transcript.expect_output(ALPS_TEST_SOURCE_DIR "/boostdump2.output");
}
