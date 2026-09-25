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
#include "../ctint/update_kernel.hpp"

void InteractionExpansion::interaction_expansion_step(void)
{ return alps::ctint::update_kernel::interaction_expansion_step(*this); }

void InteractionExpansion::reset_perturbation_series(void)
{ return alps::ctint::update_kernel::reset_perturbation_series(*this); }
