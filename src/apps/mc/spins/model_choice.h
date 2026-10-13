/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2009 by Matthias Troyer <troyer@comp-phys.org>,
*                       Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

// Which classical spin simulation the parameters select. The choice was moved
// unchanged out of SpinFactory (spinmc_factory.C), together with its helpers
// countElements, findDominantMatrixString and produceError: same reads and
// same order. produce_error now names the parameters it rejects. It depends on
// parameters only, so it can be checked without a lattice.

#ifndef ALPS_APPLICATIONS_MC_SPIN_MODEL_CHOICE_H_
#define ALPS_APPLICATIONS_MC_SPIN_MODEL_CHOICE_H_

#include <alps/parameter.h>
#include <boost/config.hpp>
#include <boost/throw_exception.hpp>
#include <algorithm>
#include <stdexcept>
#include <string>

namespace spinmc {

// Moment and matrix type of the simulation; see SpinFactory::make_worker.
enum class simulation_kind {
  ising,                 // IsingMoment, MIdMatrix<double,1>
  o4,                    // ONMoment<4>, MIdMatrix<double,4>
  heisenberg_identity,   // ONMoment<3>, MIdMatrix<double,3>
  heisenberg_diagonal,   // ONMoment<3>, DiagMatrix<double,3>
  heisenberg_square,     // ONMoment<3>, SquareMatrix<double,3>
  xy_identity,           // XYMoment, MIdMatrix<double,2>
  xy_diagonal,           // XYMoment, DiagMatrix<double,2>
  xy_square,             // XYMoment, SquareMatrix<double,2>
  potts3,                // PottsMoment<3>, MIdMatrix<double,2>
  potts4,                // PottsMoment<4>, MIdMatrix<double,2>
  potts10                // PottsMoment<10>, MIdMatrix<double,2>
};

/**
 * Returns the number of white-space separated elements in the string.
 *
 * \param str the string to be analyzed.
 */
inline int count_elements(const std::string& str)
{
  int counter = 0;
  int pos = 0;
  if (str[0] == '\0') return 0;
  while (true) {
    while (str[pos] == ' ') pos++;
    if (str[pos] == '\0') return counter;
    counter++;
    while ((str[pos] !='\0') && (str[pos] != ' ')) pos++;
  }
}

/**
 * Returns the number of elements in the string with the largest number of elements.
 */
inline int find_dominant_matrix_string(const alps::Parameters& parms)
{
  int count = 0;
  /* should be done more elegantly ... */
  if (parms.defined("J"))
    count = std::max BOOST_PREVENT_MACRO_SUBSTITUTION (count, count_elements(parms["J"]));
  if (parms.defined("J'"))
    count = std::max BOOST_PREVENT_MACRO_SUBSTITUTION (count, count_elements(parms["J'"]));
  if (parms.defined("J0"))
    count = std::max BOOST_PREVENT_MACRO_SUBSTITUTION (count, count_elements(parms["J0"]));
  if (parms.defined("J1"))
    count = std::max BOOST_PREVENT_MACRO_SUBSTITUTION (count, count_elements(parms["J1"]));
  if (parms.defined("J2"))
    count = std::max BOOST_PREVENT_MACRO_SUBSTITUTION (count, count_elements(parms["J2"]));
  if (parms.defined("J3"))
    count = std::max BOOST_PREVENT_MACRO_SUBSTITUTION (count, count_elements(parms["J3"]));
  if (parms.defined("D"))
    count = std::max BOOST_PREVENT_MACRO_SUBSTITUTION (count, count_elements(parms["D"]));
  if (parms.defined("D'"))
    count = std::max BOOST_PREVENT_MACRO_SUBSTITUTION (count, count_elements(parms["D'"]));
  if (parms.defined("D0"))
    count = std::max BOOST_PREVENT_MACRO_SUBSTITUTION (count, count_elements(parms["D0"]));
  if (parms.defined("D1"))
    count = std::max BOOST_PREVENT_MACRO_SUBSTITUTION (count, count_elements(parms["D1"]));
  if (parms.defined("D2"))
    count = std::max BOOST_PREVENT_MACRO_SUBSTITUTION (count, count_elements(parms["D2"]));
  if (parms.defined("D3"))
    count = std::max BOOST_PREVENT_MACRO_SUBSTITUTION (count, count_elements(parms["D3"]));

  return count;
}

/**
 * Reports an invalid number of elements in the input for a matrix parameter,
 * naming the parameters with the offending number of values. It reads only
 * parameters that are defined, so a missing J or D is not reported instead.
 */
inline void produce_error(const alps::Parameters& parms, int count)
{
  static const char* const counted[] = {"J", "J'", "J0", "J1", "J2", "J3",
                                        "D", "D'", "D0", "D1", "D2", "D3"};
  std::string message = "invalid parameters: cannot produce model "
    + std::string(parms["MODEL"]);
  if (count == 0)
    message += " without J, J', J0-J3, D, D' or D0-D3";
  else {
    message += " from " + std::to_string(count) + " values in";
    for (const char* name : counted)
      if (parms.defined(name) && count_elements(parms[name]) == count)
        message += std::string(" ") + name + "=\"" + std::string(parms[name]) + "\"";
  }
  boost::throw_exception(std::runtime_error(message));
}

/**
 * The type of the produced worker depends on the model given in the parameter
 * file and on the number of input elements for the matrix parameters. The
 * chosen matrix type is the most general one (over all matrix parameters) of
 * all most specialized ones that are compatible with the user input.
 */
inline simulation_kind choose_simulation(const alps::Parameters& parms)
{
  int maxElemCount = find_dominant_matrix_string(parms);
  if (parms["MODEL"]=="Ising")
    if (maxElemCount == 1)
      return simulation_kind::ising;
    else {
      produce_error(parms, maxElemCount);
    }
  else if (parms["MODEL"]=="O(4)")
    return simulation_kind::o4;
  else if (parms["MODEL"]=="Heisenberg")
    switch (maxElemCount) {
      case 1:
        return simulation_kind::heisenberg_identity;
      case 3:
        return simulation_kind::heisenberg_diagonal;
      case 6:
        return simulation_kind::heisenberg_square;
      case 9:
        return simulation_kind::heisenberg_square;
      default:
        produce_error(parms, maxElemCount);
    }
  else if (parms["MODEL"]=="XY")
    switch (maxElemCount) {
      case 1:
        return simulation_kind::xy_identity;
      case 2:
        return simulation_kind::xy_diagonal;
      case 3:
        return simulation_kind::xy_square;
      case 4:
        return simulation_kind::xy_square;
      default:
        produce_error(parms, maxElemCount);
    }

  else if (parms["MODEL"]=="Potts")
    switch (int(parms["q"])) {
    case 3:
      return simulation_kind::potts3;
    case 4:
      return simulation_kind::potts4;
    case 10:
      return simulation_kind::potts10;
    default:
      boost::throw_exception(std::runtime_error(std::string(parms["q"])
            + "-state Potts model not implemented in file factory.C"));
    }
  else
    boost::throw_exception(std::runtime_error(std::string(parms["MODEL"])
              + " model not implemented in file factory.C"));
  // not reached: every other path returns or throws
  std::abort();
}

} // namespace spinmc

#endif
