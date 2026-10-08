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
#include <alps/expression/symbol_table.h>
#include <vector>
#include <iostream>
#include <string>

#ifdef BOOST_NO_ARGUMENT_DEPENDENT_LOOKUP
using namespace alps;
#endif

int main()
{
#ifndef BOOST_NO_EXCEPTIONS
  try {
#endif

    // The cases share MODEL = "spin", J = 1 and L = 4.
    const std::vector<std::vector<std::pair<std::string, std::string> > > cases{
      {{"LATTICE", "square lattice"}},
      {{"LATTICE", "frustrated square lattice"}},
      {{"J'", "1"}, {"LATTICE", "frustrated square lattice"}},
      {{"J'", "-1"}, {"LATTICE", "frustrated square lattice"}},
      {{"LATTICE", "triangular lattice"}},
      {{"LATTICE", "chain lattice"}},
      {{"Gamma", "1"}, {"LATTICE", "chain lattice"}},
      {{"J", "-1"}, {"Gamma", "1"}, {"LATTICE", "chain lattice"}}};
    for (std::size_t i=0;i<cases.size();++i) {
      alps::SymbolTable parms;
      parms["MODEL"] = "spin";
      parms["J"] = 1;
      parms["L"] = 4;
      for (auto const& [key, value] : cases[i]) parms[key] = value;
      alps::ModelLibrary models(parms);
      alps::graph_helper<> lattice(parms);
      alps::HamiltonianDescriptor<short> ham(models.get_hamiltonian(lattice,parms,true));
      std::cout << ham;
    }

#ifndef BOOST_NO_EXCEPTIONS
}
catch (std::exception& exc) {
  std::cerr << exc.what() << "\n";
  return -1;
}
catch (...) {
  std::cerr << "Fatal Error: Unknown Exception!\n";
  return -2;
}
#endif
  return 0;
}
