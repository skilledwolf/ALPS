/****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2012 by Hartmut Hafermann <hafermann@cpht.polytechnique.fr>
 *                       Emanuel Gull <gull@pks.mpg.de>,
 *
 *  based on an earlier version by Philipp Werner and Emanuel Gull
 *
 *
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
 *
 *****************************************************************************/
#ifndef HYB_EVALUATE
#define HYB_EVALUATE
#include "hyb.hpp"
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/bessel.hpp>

void evaluate_basics(const hybridization::results_type &results, const alps::params &parms, const alps::params &output, alps::hdf5::archive &solver_output);
void evaluate_time(const hybridization::results_type &results, const alps::params &parms, const alps::params &output, alps::hdf5::archive &solver_output);
void evaluate_freq(const hybridization::results_type &results, const alps::params &parms, const alps::params &output, alps::hdf5::archive &solver_output);
void evaluate_legendre(const hybridization::results_type &results, const alps::params &parms, const alps::params &output, alps::hdf5::archive &solver_output);
void evaluate_nnt(const hybridization::results_type &results, const alps::params &parms, const alps::params &output, alps::hdf5::archive &solver_output);
void evaluate_nnw(const hybridization::results_type &results, const alps::params &parms, const alps::params &output, alps::hdf5::archive &solver_output);
void evaluate_sector_statistics(const hybridization::results_type &results, const alps::params &parms, const alps::params &output, alps::hdf5::archive &solver_output);
void evaluate_2p(const hybridization::results_type &results, const alps::params &parms, const alps::params &output, alps::hdf5::archive &solver_output);

inline std::complex<double> t(int n, int l){//transformation matrix from Legendre to Matsubara basis
  std::complex<double> i_c(0., 1.);
  return (std::sqrt(static_cast<double>(2*l+1))/std::sqrt(static_cast<double>(2*n+1))) * std::exp(i_c*(n+0.5)*boost::math::constants::pi<double>()) * std::pow(i_c,l) * boost::math::cyl_bessel_j(l+0.5,(n+0.5)*boost::math::constants::pi<double>());
}

#endif
