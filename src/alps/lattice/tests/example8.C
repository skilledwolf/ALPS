#include <alps/testing/stream_fixture.hpp>
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

/* $Id$ */

#include <alps/lattice.h>
#include <iostream>
#include <fstream>

#ifdef BOOST_NO_ARGUMENT_DEPENDENT_LOOKUP
using namespace alps;
#endif

TEST(LatticeSerialization, Example8) {
    alps::testing::StreamFixture transcript(ALPS_TEST_SOURCE_DIR "/example8.input");
    { // Flush serialization objects before checking the captured stream.


#ifndef BOOST_NO_EXCEPTIONS
  try {
#endif

    // read parameters
    alps::Parameters parameters;
    std::cin >> parameters;
    // create a graph factory with default graph type
    alps::graph_helper<> lattice(parameters);
    // write the graph created from the input in XML
    std::cout << lattice.graph();

#ifndef BOOST_NO_EXCEPTIONS
}
catch (std::exception& e)
{
  std::cerr << "Caught exception: " << e.what() << "\n";
  FAIL() << "Unexpected exception in serialization contract";
}
catch (...)
{
  std::cerr << "Caught unknown exception\n";
  FAIL() << "Unexpected exception in serialization contract";
}
#endif
    }
  transcript.expect_output(ALPS_TEST_SOURCE_DIR "/example8.output");
}
