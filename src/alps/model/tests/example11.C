#include <alps/testing/stream_fixture.hpp>
/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2003-2004 by Matthias Troyer <troyer@itp.phys.ethz.ch>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include <alps/model.h>
#include <iostream>

TEST(ModelSerialization, Example11) {
    alps::testing::StreamFixture transcript(ALPS_TEST_SOURCE_DIR "/example11.input");
    { // Flush serialization objects before checking the captured stream.


#ifndef BOOST_NO_EXCEPTIONS
  try {
#endif
    alps::Parameters parms;
    std::cin >> parms;
    alps::ModelLibrary models(parms);
    alps::graph_helper<> lattices(parms);
    alps::HamiltonianDescriptor<short> ham(models.get_hamiltonian(parms["MODEL"]));
    parms.copy_undefined(ham.default_parameters());
    ham.set_parameters(parms);
    alps::basis_states_descriptor<short> basis(ham.basis(),lattices.graph());
    for (int i=0;i<basis.get_basis().constraints().size();++i)
      std::cout << "Constraint: " << basis.get_basis().constraints()[i].first << "=" 
                <<  basis.get_basis().constraints()[i].second << std::endl;

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
  transcript.expect_output(ALPS_TEST_SOURCE_DIR "/example11.output");
}
