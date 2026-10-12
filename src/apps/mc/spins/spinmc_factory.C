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
#include "model_choice.h"

#include "spinrun.h"

#include "ising.h"
#include "on.h"
#include "potts.h"

void SpinFactory::print_copyright(std::ostream& out) const
{
  out << "Generic classical Monte Carlo program using local or cluster updates\n"
      << "  available from http://alps.comp-phys.org/\n"
      << "  copyright(c) 1999-2007 by Matthias Troyer <troyer@comp-phys.org>\n"
      << "                            Mathias Koerner <mkoerner@comp-phys.org>\n"
      << "\n";
}


alps::scheduler::MCSimulation* SpinFactory::make_task(
        const alps::ProcessList& w,
        const boost::filesystem::path& fn) const
{
  return new alps::scheduler::MCSimulation(w,fn);
}

/**
 * Creates a worker, according to the parameters. spinmc::choose_simulation
 * (model_choice.h) selects the moment and matrix types.
 */
alps::scheduler::Worker* SpinFactory::make_worker(
        const alps::ProcessList& where,
        const alps::Parameters& parms,int node) const
{
  switch (spinmc::choose_simulation(parms)) {
    case spinmc::simulation_kind::ising:
      return new SpinRun<IsingMoment,MIdMatrix<double,1> >(where,parms,node);
    case spinmc::simulation_kind::o4:
      return new SpinRun<ONMoment<4>,MIdMatrix<double,4> >(where,parms,node);
    case spinmc::simulation_kind::heisenberg_identity:
      return new SpinRun< ONMoment<3>, MIdMatrix<double,3> >(where, parms, node);
    case spinmc::simulation_kind::heisenberg_diagonal:
      return new SpinRun< ONMoment<3>, DiagMatrix<double,3> >(where, parms, node);
    case spinmc::simulation_kind::heisenberg_square:
      return new SpinRun< ONMoment<3>, SquareMatrix<double,3> >(where, parms, node);
    case spinmc::simulation_kind::xy_identity:
      return new SpinRun< XYMoment, MIdMatrix<double,2> >(where, parms, node);
    case spinmc::simulation_kind::xy_diagonal:
      return new SpinRun< XYMoment, DiagMatrix<double,2> >(where, parms, node);
    case spinmc::simulation_kind::xy_square:
      return new SpinRun< XYMoment, SquareMatrix<double,2> >(where, parms, node);
    case spinmc::simulation_kind::potts3:
      return new SpinRun< PottsMoment<3>, MIdMatrix<double,2> >(where,parms,node);
    case spinmc::simulation_kind::potts4:
      return new SpinRun< PottsMoment<4>, MIdMatrix<double,2> >(where,parms,node);
    case spinmc::simulation_kind::potts10:
      return new SpinRun< PottsMoment<10>, MIdMatrix<double,2> >(where,parms,node);
  }
  return 0;
}
