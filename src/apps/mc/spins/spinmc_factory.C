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

/* $Id$ */

#include "factory.h"
#include "abstractspinsim.h"

#include "spinsim.h"

#include "ising.h"
#include "on.h"
#include "potts.h"

namespace {
template<class Moment>
alps::scheduler::Worker* make_vector_spin(int count, alps::ProcessList const& where,
                                         alps::Parameters const& parms, int node) {
  constexpr int dim = Moment::dim;
  switch (count) {
    case 1: return new SpinSim<Moment, MIdMatrix<double, dim>>(where, parms, node);
    case dim: return new SpinSim<Moment, DiagMatrix<double, dim>>(where, parms, node);
    case dim * (dim + 1) / 2:
    case dim * dim: return new SpinSim<Moment, SquareMatrix<double, dim>>(where, parms, node);
    default: return nullptr;
  }
}
}

void SpinFactory::print_copyright(std::ostream& out) const
{
  out << "Generic classical Monte Carlo program using local or cluster updates\n"
      << "  available from http://alps.comp-phys.org/\n"
      << "  copyright(c) 1999-2007 by Matthias Troyer <troyer@comp-phys.org>\n"
      << "                            Mathias Koerner <mkoerner@comp-phys.org>\n"
      << " for details see the publication:\n"
      << " A.F. Albuquerque et al., J. of Magn. and Magn. Materials 310, 1187 (2007).\n\n";
}


alps::scheduler::MCSimulation* SpinFactory::make_task(
        const alps::ProcessList& w,
        const boost::filesystem::path& fn) const
{
  return new alps::scheduler::MCSimulation(w,fn);
}

/** 
 * Returns the number of white-space separated elements in the string.
 *
 * \param str the string to be analyzed.
 */
int SpinFactory::countElements(const std::string& str) const
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
int SpinFactory::findDominantMatrixString(const alps::Parameters& parms) const
{
  int count = 0;
  for (auto key : {"J", "J'", "J0", "J1", "J2", "J3", "D", "D'", "D0", "D1", "D2", "D3"})
    if (parms.defined(key)) count = (std::max)(count, countElements(parms[key]));

  return count;
}

/**
 * Outputs an error message and aborts the computation if an invalid number of
 * elements is provided in the input for a matrix parameter.
 */
void SpinFactory::produceError(const alps::Parameters& parms) const
{
  std::cerr << "Invalid parameters given:\n"
            << "Cannot produce model " << parms["MODEL"] << " with strings "
            << parms["D"] << " for D and " << parms["J"] 
            << " for J.\nAborting\n";
  boost::throw_exception(std::runtime_error("invalid parameters"));
  std::exit(-5);
}

/**
 * Creates a worker, according to the parameters. The type of the produced 
 * worker depends on the model given in the parameter file and on the number of
 * input elements for the matrix parameters. The chosen matrix type is the most
 * general one (over all matrix parameters) of all most specialized ones that
 *  are compatible with the user input.
 */
alps::scheduler::Worker* SpinFactory::make_worker(
        const alps::ProcessList& where,
        const alps::Parameters& parms,int node) const
{
  int maxElemCount = findDominantMatrixString(parms);
  if (parms["MODEL"]=="Ising") {
    if (maxElemCount == 1)
      return new SpinSim<IsingMoment,MIdMatrix<double,1> >(where,parms,node);
  } else if (parms["MODEL"]=="O(4)") {
    return new SpinSim<ONMoment<4>,MIdMatrix<double,4> >(where,parms,node);
  } else if (parms["MODEL"]=="Heisenberg") {
    if (auto* worker = make_vector_spin<ONMoment<3>>(maxElemCount, where, parms, node)) return worker;
  } else if (parms["MODEL"]=="XY") {
    if (auto* worker = make_vector_spin<XYMoment>(maxElemCount, where, parms, node)) return worker;
  }
  else if (parms["MODEL"]=="Potts")
    switch (int(parms["q"])) {
    case 3:
      return new SpinSim< PottsMoment<3>, MIdMatrix<double,2> >(where,parms,node);
    case 4:
      return new SpinSim< PottsMoment<4>, MIdMatrix<double,2> >(where,parms,node);
    case 10:
      return new SpinSim< PottsMoment<10>, MIdMatrix<double,2> >(where,parms,node);
    default:
      boost::throw_exception(std::runtime_error(std::string(parms["q"]) 
            + "-state Potts model not implemented in file factory.C"));
    }
  else 
    boost::throw_exception(std::runtime_error(std::string(parms["MODEL"])
              + " model not implemented in file factory.C"));
  produceError(parms);
  return 0;
}
