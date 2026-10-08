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
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#ifdef BOOST_NO_ARGUMENT_DEPENDENT_LOOKUP
using namespace alps;
#endif

int main() {
#ifndef BOOST_NO_EXCEPTIONS
  try {
#endif
    typedef alps::coordinate_graph_type graph_t;
    typedef alps::parity_t parity_t;

    // each lattice, optionally with its BACKBONE_TYPES
    const std::vector<std::pair<std::string, char const*> > cases{
      {"square lattice", nullptr}, {"frustrated square lattice", nullptr},
      {"frustrated square lattice", ""}, {"frustrated square lattice", "0"},
      {"frustrated square lattice", "1"}, {"frustrated square lattice", "0 1"}};
    for (auto const& [name, backbone] : cases) {
      alps::params p;
      p["L"] = 4;
      p["LATTICE"] = name;
      if (backbone) p["BACKBONE_TYPES"] = std::string(backbone);
      // create the lattice

      alps::graph_helper<> lattice(p);
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
    exit(-1);
  }
  catch (...) {
    std::cerr << "Caught unknown exception\n";
    exit(-2);
  }
#endif
  return 0;
}
