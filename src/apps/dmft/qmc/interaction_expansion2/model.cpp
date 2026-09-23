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

double HalfFillingHubbardInteractionExpansion::try_add()
{ return alps::ctint::update_kernel::half_filling_try_add(*this); }

void HalfFillingHubbardInteractionExpansion::perform_add()
{ return alps::ctint::update_kernel::half_filling_perform_add(*this); }

void HalfFillingHubbardInteractionExpansion::reject_add()
{ return alps::ctint::update_kernel::half_filling_reject_add(*this); }

double HalfFillingHubbardInteractionExpansion::try_remove(unsigned int vertex_nr)
{ return alps::ctint::update_kernel::half_filling_try_remove(*this, vertex_nr); }

void HalfFillingHubbardInteractionExpansion::perform_remove(unsigned int vertex_nr)
{ return alps::ctint::update_kernel::half_filling_perform_remove(*this, vertex_nr); }

void HalfFillingHubbardInteractionExpansion::reject_remove()
{ return alps::ctint::update_kernel::half_filling_reject_remove(*this); }

double HubbardInteractionExpansion::try_add()
{ return alps::ctint::update_kernel::hubbard_try_add(*this); }

void HubbardInteractionExpansion::perform_add()
{ return alps::ctint::update_kernel::hubbard_perform_add(*this); }

void HubbardInteractionExpansion::reject_add()
{ return alps::ctint::update_kernel::hubbard_reject_add(*this); }

double HubbardInteractionExpansion::try_remove(unsigned int vertex_nr)
{ return alps::ctint::update_kernel::hubbard_try_remove(*this, vertex_nr); }

void HubbardInteractionExpansion::perform_remove(unsigned int vertex_nr)
{ return alps::ctint::update_kernel::hubbard_perform_remove(*this, vertex_nr); }

void HubbardInteractionExpansion::reject_remove()
{ return alps::ctint::update_kernel::hubbard_reject_remove(*this); }

double MultiBandDensityHubbardInteractionExpansion::try_add()
{ return alps::ctint::update_kernel::multiband_try_add(*this); }

void MultiBandDensityHubbardInteractionExpansion::perform_add()
{ return alps::ctint::update_kernel::multiband_perform_add(*this); }

void MultiBandDensityHubbardInteractionExpansion::reject_add()
{ return alps::ctint::update_kernel::multiband_reject_add(*this); }

double MultiBandDensityHubbardInteractionExpansion::try_remove(unsigned int vertex_nr)
{ return alps::ctint::update_kernel::multiband_try_remove(*this, vertex_nr); }

void MultiBandDensityHubbardInteractionExpansion::perform_remove(unsigned int vertex_nr)
{ return alps::ctint::update_kernel::multiband_perform_remove(*this, vertex_nr); }

void MultiBandDensityHubbardInteractionExpansion::reject_remove()
{ return alps::ctint::update_kernel::multiband_reject_remove(*this); }
