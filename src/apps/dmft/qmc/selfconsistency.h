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

/* $Id: selfconsistency.h 370 2009-08-05 10:08:34Z fuchs $ */

#ifndef ALPS_DMFT_SELFCONSISTENCY_H
#define ALPS_DMFT_SELFCONSISTENCY_H

/// @file selfconsistency.h
/// @brief declares the selfconsistency loop functions

#include "solver.h"
#include <alps/run_config.hpp>
#include "hilberttransformer.h"
#include "fouriertransform.h"
#include "green_function.h"

/// performs a DMFT self-consistency loop in imaginary time for a solver that
/// receives the Bethe-lattice hybridization function
///
/// @param solver   The impurity solver. It is left in the state after the final iteration to 
///                 retrieve additional information.
/// @param hilbert  symmetrizes the dressed Green function.
/// @param initial  the initial Green function.
extern void F_selfconsistency_loop(alps::run_configuration& run, ImpuritySolver& solver, HilbertTransformer& hilbert,
                                  itime_green_function_t initial);


/// performs a DMFT self-consistency loop until convergence
///
/// @param parms The input parameters for the simulation.
/// @param solver   An impurity solver that takes both the bare GF in imaginary time AND in Matsubara freqencies, but returns
/// @param hilbert  a Hilbert transformation object. It performs the Hilbert transformation, taking its arguments in Frequency space.
void selfconsistency_loop_omega(alps::run_configuration& run, MatsubaraImpuritySolver& solver,
                               FrequencySpaceHilbertTransformer& hilbert, matsubara_green_function_t initial);


//void selfconsistency_loop_DCA(const alps::run_configuration& run, MatsubaraImpuritySolver& solver, DCATransformer& clustertrans);

#endif
