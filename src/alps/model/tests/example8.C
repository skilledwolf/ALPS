#include <alps/testing/stream_fixture.hpp>
/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2003-2006 by Matthias Troyer <troyer@itp.phys.ethz.ch>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include <alps/model.h>
#include <alps/lattice.h>
#include <alps/parameter.h>
#include <iostream>
#include <string>

#ifdef BOOST_NO_ARGUMENT_DEPENDENT_LOOKUP
using namespace alps;
#endif

TEST(ModelSerialization, Example8) {
    alps::testing::StreamFixture transcript(ALPS_TEST_SOURCE_DIR "/example8.input");
    { // Flush serialization objects before checking the captured stream.

#ifndef BOOST_NO_EXCEPTIONS
  try {
#endif

    alps::ParameterList parms;
    std::cin >> parms;
    for (int i=0;i<parms.size();++i) {
      alps::ModelLibrary models(parms[i]);
      alps::graph_helper<> lattice(parms[i]);
      alps::HamiltonianDescriptor<short> ham(models.get_hamiltonian(lattice,parms[i]));
      parms[i].copy_undefined(ham.default_parameters());
      ham.set_parameters(parms[i]);
      if (has_sign_problem(ham,lattice,parms[i]))
        std::cout << "Model " << i+1 << " has a sign problem.\n";
      else
        std::cout << "Model " << i+1 << " has no sign problem.\n";
    }

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
  transcript.expect_output(ALPS_TEST_SOURCE_DIR "/example8.output");
}
