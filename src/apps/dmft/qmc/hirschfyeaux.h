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

/* $Id: hirschfyeaux.h 328 2008-09-09 19:02:18Z gullc $ */

/// @file hirschfyesim.h
/// @brief the actual Hirsch-Fye simulation

#ifndef ALPS_DMFT_HIRSCHFYEAUX_H
#define ALPS_DMFT_HIRSCHFYEAUX_H

#include <alps/config.h> // needed to set up correct bindings
#include <boost/numeric/bindings/ublas.hpp>
#include <boost/numeric/bindings/blas.hpp>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>

typedef boost::numeric::ublas::matrix<double,boost::numeric::ublas::column_major> dense_matrix;

// acosh(exp(x)), evaluated without squaring exp(x) or cancellation at x=0.
// The rank-one update also requires exp(2*lambda) to be representable.
inline double hirschfye_lambda(double beta, double u, std::uint64_t slices) {
  if (!(beta > 0.) || !(u >= 0.) || !slices || !std::isfinite(beta) || !std::isfinite(u))
    throw std::invalid_argument("Hirsch-Fye Hubbard-Stratonovich parameters are invalid");
  const double x = (beta / slices) * (u / 2.);
  const double lambda = x + std::log1p(std::sqrt(-std::expm1(-2. * x)));
  if (!std::isfinite(lambda) || !std::isfinite(std::exp(2. * lambda)))
    throw std::invalid_argument("Hirsch-Fye Hubbard-Stratonovich exponent is not representable");
  return lambda;
}

double norm_max(dense_matrix const& matrix);
// Returns the determinant sign of the solved configuration.
int update_from_zero(dense_matrix& green, dense_matrix const& green0,
                     std::vector<int> const& spins, double lambda);


/// attempt single spin flip and (if successful) compute Green's function matrix for the new configuration
template<class RNG> 
void update_single_spin(RNG & rng, dense_matrix & Green_up, dense_matrix & Green_down, std::vector<int> & spins, double lambda, int &sign) 
{
  fortran_int_t N(spins.size());
  
  // choose random site
  int site = (unsigned int)(N*rng());
  // calculate ratio of determinants
  double r_up = 1 + (1-Green_up(site,site))*(exp(-2*lambda*spins[site])-1);
  double r_down = 1 + (1-Green_down(site,site))*(exp(2*lambda*spins[site])-1);
  double det_rat = r_up*r_down;
  if (!std::isfinite(det_rat))
    throw std::runtime_error("Hirsch-Fye determinant ratio is not finite");
  const double absolute_ratio = std::abs(det_rat);
  const double probability = absolute_ratio / (1. + absolute_ratio);
  
  // if update successful ...
  if (rng() < probability) {
    if(det_rat<0)
      sign *=-1;
    // update Green's function
    double tmp1 = exp(-2*lambda*spins[site])-1;
    double tmp2 = exp(+2*lambda*spins[site])-1;
    
    std::vector<double> vi_up(N), uj_up(N), vi_dn(N), uj_dn(N);
    double alpha_up=tmp1/(1+(1-Green_up(site,site))*tmp1);
    double alpha_dn=tmp2/(1+(1-Green_down(site,site))*tmp2);
    for(int i=0;i<N;++i){
      vi_up[i]=Green_up(i,site)-(i==(int)site);
      vi_dn[i]=Green_down(i,site)-(i==(int)site);
      uj_up[i]=Green_up  (site, i);
      uj_dn[i]=Green_down(site, i);
    }
    fortran_int_t one=1;
    FORTRAN_ID(dger)(&N,&N,&alpha_up,vi_up.data(),&one,uj_up.data(),&one,&(Green_up  (0,0)),&N);
    FORTRAN_ID(dger)(&N,&N,&alpha_dn,vi_dn.data(),&one,uj_dn.data(),&one,&(Green_down(0,0)),&N);
    // update spin
    spins[site] = -spins[site];
  }
}

#endif
