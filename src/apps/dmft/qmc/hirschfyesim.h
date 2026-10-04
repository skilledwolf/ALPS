/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2009 by Emanuel Gull <gull@phys.columbia.edu>
 *                              Philipp Werner <werner@itp.phys.ethz.ch>,
 *                              Matthias Troyer <troyer@comp-phys.org>
 *
 *
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
 *
 *****************************************************************************/

#ifndef ALPS_DMFT_HIRSCHFYESIM_H
#define ALPS_DMFT_HIRSCHFYESIM_H

#include "hirschfyeaux.h"
#include "green_function.h"
#include <alps/run_config.hpp>
#include <alps/ngs/random01.hpp>
#include <alps/alea/batch.hpp>
#include <array>
#include <map>
#include <string>

// Private solver engine. The standalone driver owns stopping and publication;
// these measurement snapshots do not contain a resumable solver state.
class HirschFyeRun {
public:
  using results_type = std::map<std::string, alps::alea::batch_result<double>>;
  HirschFyeRun(alps::run_configuration const&, matsubara_green_function_t const& g0, int rank = 0);
  void dostep();
  double work_done() const;
  results_type collect_results(alps::alea::reducer const* reduction = nullptr) const;
  static std::pair<matsubara_green_function_t, itime_green_function_t>
    get_result(results_type const&, alps::params const&);

  void green_vector_from_matrix(itime_green_function_t&, dense_matrix const&, dense_matrix const&) const;
  void green_matrix_from_vector(itime_green_function_t const&, dense_matrix&, dense_matrix&) const;

protected:
  // Unsigned Sign, then each spin's signed physical G(tau)[0..N] + sign.
  std::array<alps::alea::batch_acc<double>, 3> measurements;
  alps::params parameters;
  alps::random01 random_01;
  std::uint64_t sweeps = 0, thermalization_sweeps, total_sweeps;
  double beta;
  int N;
  double lambda;
  std::uint64_t N_check = 1000, check_counter = 0;
  double tolerance;
  dense_matrix Green0_up, Green_up, Green0_down, Green_down;
  std::vector<int> spins;
  int sign = 1;
  itime_green_function_t bare_green_tau, green_tau;
};

#endif
