/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2009 by Emanuel Gull <gull@phys.columbia.edu>
 *                              Philipp Werner <werner@itp.phys.ethz.ch>,
 *                              Sebastian Fuchs <fuchs@theorie.physik.uni-goettingen.de>
 *                              Matthias Troyer <troyer@comp-phys.org>
 *
 *
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include "interaction_expansion.hpp"

void InteractionExpansionRun::load(alps::IDump &/*dump*/)
{
  std::cerr<<"loading: doing nothing."<<std::endl;
}

void InteractionExpansionRun::save(alps::ODump &/*dump*/) const
{
  std::cerr<<"saving: doing nothing"<<std::endl;
}
