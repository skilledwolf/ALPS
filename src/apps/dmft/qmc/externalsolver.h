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

/* $Id: externalsolver.h 157 2006-04-04 15:11:22Z gullc $ */

#ifndef ALPS_DMFT_EXTERNALSOLVER_H
#define ALPS_DMFT_EXTERNALSOLVER_H

/// @file externalsolver.h
/// @brief declares the external solver
/// @sa ExternalSolver

#include "solver.h"
#include "run_config.h"
#include <filesystem>

/// Calls an impurity solver with one TOML run file. Scientific arrays travel
/// separately in HDF5 as /Delta_<flavor> or /G0_<flavor> vectors; the solver
/// writes /G_tau and optionally /G_omega in the DMFT Green-function layout.
class ExternalSolver 
 : public ImpuritySolver
 , public MatsubaraImpuritySolver 
{
public:
    explicit ExternalSolver(const alps::run_configuration& configuration);
  
    ImpuritySolver::result_type solve(
              const itime_green_function_t& G0
            , const alps::params& parms);
    
    MatsubaraImpuritySolver::result_type solve_omega(
              const matsubara_green_function_t& G0_omega
            , const alps::params& parms);
private:
    struct invocation_files;
    alps::run_configuration solver_run(const alps::params& parameters, const invocation_files& files) const;
    void write_hybridization(alps::run_configuration& run, const alps::params& parameters,
                             const itime_green_function_t& delta, alps::hdf5::archive& input) const;
    void call(alps::run_configuration run, const invocation_files& files) const;

    alps::run_configuration configuration_;
    alps::dmft::solver_kind kind_;
    std::string schema_;
    bool delta_;
    std::filesystem::path executable_;
};



#endif
