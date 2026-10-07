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

#include <alps/parser/xmlparser.h>
#include <iostream>
#include <cstdlib>

TEST(XmlSerialization, Xmlparser) {
    alps::testing::StreamFixture transcript(ALPS_TEST_SOURCE_DIR "/xmlparser.input");
    { // Flush serialization objects before checking the captured stream.

#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  alps::PrintXMLHandler handler;
  alps::XMLParser parser(handler);
  
  parser.parse(std::cin);

#ifndef BOOST_NO_EXCEPTIONS
}
catch (std::exception& exp) {
  std::cerr << exp.what() << std::endl;
  FAIL() << "Unexpected exception in serialization contract";
}
#endif
    }
  transcript.expect_output(ALPS_TEST_SOURCE_DIR "/xmlparser.output");
}
