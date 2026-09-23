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

double HalfFillingHubbardInteractionExpansionRun::try_add()
{ return alps::ctint::update_kernel::half_filling_try_add(*this); }

void HalfFillingHubbardInteractionExpansionRun::perform_add()
{ return alps::ctint::update_kernel::half_filling_perform_add(*this); }

void HalfFillingHubbardInteractionExpansionRun::reject_add()
{ return alps::ctint::update_kernel::half_filling_reject_add(*this); }

double HalfFillingHubbardInteractionExpansionRun::try_remove(unsigned int vertex_nr)
{ return alps::ctint::update_kernel::half_filling_try_remove(*this, vertex_nr); }

void HalfFillingHubbardInteractionExpansionRun::perform_remove(unsigned int vertex_nr)
{ return alps::ctint::update_kernel::half_filling_perform_remove(*this, vertex_nr); }

void HalfFillingHubbardInteractionExpansionRun::reject_remove()
{ return alps::ctint::update_kernel::half_filling_reject_remove(*this); }

double HubbardInteractionExpansionRun::try_add()
{ return alps::ctint::update_kernel::hubbard_try_add(*this); }

void HubbardInteractionExpansionRun::perform_add()
{ return alps::ctint::update_kernel::hubbard_perform_add(*this); }

void HubbardInteractionExpansionRun::reject_add()
{ return alps::ctint::update_kernel::hubbard_reject_add(*this); }

double HubbardInteractionExpansionRun::try_remove(unsigned int vertex_nr)
{ return alps::ctint::update_kernel::hubbard_try_remove(*this, vertex_nr); }

void HubbardInteractionExpansionRun::perform_remove(unsigned int vertex_nr)
{ return alps::ctint::update_kernel::hubbard_perform_remove(*this, vertex_nr); }

void HubbardInteractionExpansionRun::reject_remove()
{ return alps::ctint::update_kernel::hubbard_reject_remove(*this); }

double MultiBandDensityHubbardInteractionExpansionRun::try_add()
{ return alps::ctint::update_kernel::multiband_try_add(*this); }

void MultiBandDensityHubbardInteractionExpansionRun::perform_add()
{ return alps::ctint::update_kernel::multiband_perform_add(*this); }

void MultiBandDensityHubbardInteractionExpansionRun::reject_add()
{ return alps::ctint::update_kernel::multiband_reject_add(*this); }

double MultiBandDensityHubbardInteractionExpansionRun::try_remove(unsigned int vertex_nr)
{ return alps::ctint::update_kernel::multiband_try_remove(*this, vertex_nr); }

void MultiBandDensityHubbardInteractionExpansionRun::perform_remove(unsigned int vertex_nr)
{ return alps::ctint::update_kernel::multiband_perform_remove(*this, vertex_nr); }

void MultiBandDensityHubbardInteractionExpansionRun::reject_remove()
{ return alps::ctint::update_kernel::multiband_reject_remove(*this); }
