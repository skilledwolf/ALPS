#include <alps/testing/stream_fixture.hpp>
/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2003 by Matthias Troyer <troyer@itp.phys.ethz.ch>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include <alps/parser/xslt_path.h>
#include <alps/model.h>
#include <fstream>
#include <iostream>

TEST(ModelSerialization, Example1) {
    alps::testing::StreamFixture transcript;
    { // Flush serialization objects before checking the captured stream.


#ifndef BOOST_NO_EXCEPTIONS
  try {
#endif
    // create the library from an XML file
    std::ifstream in(alps::search_xml_library_path("models.xml"));
    alps::ModelLibrary lib(in);

    // write the library in XML
    std::cout << lib;
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
  transcript.expect_output(ALPS_TEST_SOURCE_DIR "/example1.output");
}
