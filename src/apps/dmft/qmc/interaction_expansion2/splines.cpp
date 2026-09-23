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

double InteractionExpansion::green0_spline(const creator &cdagger, const annihilator &c) const
{ return alps::ctint::update_kernel::green0_spline(*this, cdagger, c); }

double InteractionExpansion::green0_spline(const itime_t delta_t, const spin_t flavor, const site_t site1, const site_t site2) const
{ return alps::ctint::update_kernel::green0_spline(*this, delta_t, flavor, site1, site2); }

double InteractionExpansion::green0_spline(const itime_t delta_t, const spin_t flavor) const
{ return alps::ctint::update_kernel::green0_spline(*this, delta_t, flavor); }
