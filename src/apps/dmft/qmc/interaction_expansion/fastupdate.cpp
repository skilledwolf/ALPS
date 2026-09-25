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

double InteractionExpansionRun::fastupdate_up(const int flavor, bool compute_only_weight)
{ return alps::ctint::update_kernel::fastupdate_up(*this, flavor, compute_only_weight); }

double InteractionExpansionRun::fastupdate_down(const int operator_nr, const int flavor, bool compute_only_weight)
{ return alps::ctint::update_kernel::fastupdate_down(*this, operator_nr, flavor, compute_only_weight); }
