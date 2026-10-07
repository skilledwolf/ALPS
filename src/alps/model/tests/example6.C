#include <alps/testing/stream_fixture.hpp>
/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2003-2004 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Axel Grzesik <axel@th.physik.uni-bonn.de>
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
#include <string>

void write_set(const std::string& name, const alps::ModelLibrary& lib, 
               const alps::Parameters& p=alps::Parameters())
{
  alps::SiteBasisDescriptor<short> sitebasis=lib.get_site_basis(name);
  sitebasis.set_parameters(p);
  std::cout << "States of basis " << name << "=" << alps::site_basis<short>(sitebasis);
}

TEST(ModelSerialization, Example6) {
    alps::testing::StreamFixture transcript;
    { // Flush serialization objects before checking the captured stream.

#ifndef BOOST_NO_EXCEPTIONS
try {
#endif
  
  std::ifstream in(alps::search_xml_library_path("models.xml"));
  alps::ModelLibrary lib(in);

  alps::Parameters p;
  write_set("spinful boson",lib);
  p["boson_spin"]=2;
  p["NMax"]=2;
  write_set("spinful boson",lib,p);
  write_set("t-J",lib);
  write_set("alternative t-J",lib);
#ifndef BOOST_NO_EXCEPTIONS
}
catch (std::exception& exc) {
  std::cerr << exc.what() << "\n";
  FAIL() << "Unexpected exception in serialization contract";
}
catch (...) {
  std::cerr << "Fatal Error: Unknown Exception!\n";
  FAIL() << "Unexpected exception in serialization contract";
}
#endif
    }
  transcript.expect_output(ALPS_TEST_SOURCE_DIR "/example6.output");
}
