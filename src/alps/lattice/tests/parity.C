#include <alps/testing/stream_fixture.hpp>
/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2006 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include <alps/lattice.h>
#include <alps/parameter.h>
#include <fstream>
#include <iostream>

#ifdef BOOST_NO_ARGUMENT_DEPENDENT_LOOKUP
using namespace alps;
#endif

TEST(LatticeSerialization, Parity) {
    alps::testing::StreamFixture transcript(ALPS_TEST_SOURCE_DIR "/parity.input");
    { // Flush serialization objects before checking the captured stream.

#ifndef BOOST_NO_EXCEPTIONS
  try {
#endif
    typedef alps::coordinate_graph_type graph_t;
    typedef alps::parity_t parity_t;

    // read parameters
    alps::ParameterList parms;
    std::cin >> parms;
    
    for (alps::ParameterList::const_iterator p = parms.begin();
         p != parms.end(); ++p) {
      // create the lattice

      alps::graph_helper<> lattice(*p);
      const graph_t& graph = lattice.graph();
      
      std::cout << graph;
      
      for (graph_t::vertex_iterator vi = boost::vertices(graph).first;
           vi != boost::vertices(graph).second; ++vi) {
        std::cout << "vertex " << *vi << "'s parity is ";
        if (boost::get(parity_t(), graph, *vi)
            == alps::parity_traits<parity_t, graph_t>::white) {
          std::cout << "white\n";
        } else if (boost::get(alps::parity_t(), graph, *vi)
                   == alps::parity_traits<parity_t, graph_t>::black) {
          std::cout << "black\n";
        } else {
          std::cout << "undefined\n";
        }        
      }
    }
         
#ifndef BOOST_NO_EXCEPTIONS
  }
  catch (std::exception& e) {
    std::cerr << "Caught exception: " << e.what() << "\n";
    FAIL() << "Unexpected exception in serialization contract";
  }
  catch (...) {
    std::cerr << "Caught unknown exception\n";
    FAIL() << "Unexpected exception in serialization contract";
  }
#endif
    }
  transcript.expect_output(ALPS_TEST_SOURCE_DIR "/parity.output");
}
