/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2009 by Emanuel Gull <gull@phys.columbia.edu>
 *                              Philipp Werner <werner@itp.phys.ethz.ch>,
 *                              Sebastian Fuchs <fuchs@theorie.physik.uni-goettingen.de>
 *                              Matthias Troyer <troyer@comp-phys.org>
 *               2012 - 2013 by Jakub Imriska <jimriska@phys.ethz.ch>
 *
 *
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id: hilberttransformer.h 367 2009-08-05 10:04:31Z fuchs $ */

#ifndef ALPS_DMFT_HILBERTTRANSFORMER_H
#define ALPS_DMFT_HILBERTTRANSFORMER_H


/// @file hilberttransformer.h
/// @brief Hilbert transformations
///

/// declares the abstract base class and concrete realizations for the Hilbert transformations
/// @sa HilbertTransformer, SemicircleHilbertTransformer, FrequencySpaceHilbertTransformer
///

#include "types.h"
#include "bandstructure.h"
#include "fouriertransform.h"

/// @brief performs a Hilbert transformation
///
/// The imaginary-time loop solves with a hybridization solver, which applies the
/// self-consistency itself; the transformer supplies the initial G0 and symmetrizes.
class HilbertTransformer
{
public:
  itime_green_function_t symmetrize(const itime_green_function_t& G_tau, const bool symmetrization) const;
  virtual itime_green_function_t initial_G0(const alps::params& parms, const alps::params& input) const=0;
  virtual ~HilbertTransformer() {}
};



/// The imaginary-time transformer for a semicircle density of states
class SemicircleHilbertTransformer : public HilbertTransformer 
{
public:
  /// the constructor accepts the bandwidth
  explicit SemicircleHilbertTransformer(alps::params& parms)
    : bethe_parms(parms,true)
  {
    bethe_parms.set_parms(parms);
  }
  
  itime_green_function_t initial_G0(const alps::params& parms, const alps::params& input) const;
  
private:
  SemicircleBandstructure bethe_parms;
};






/// @brief performs a Hilbert transformation
///
/// The FrequencySpaceHilbertTransformer performs a Hilbert transformation for the self energy and density of states.
/// Arguments are expected to be in Matsubara Frequencies.
class FrequencySpaceHilbertTransformer{
public:
  
  virtual ~FrequencySpaceHilbertTransformer() {}
  
  /// the function call operator performs a Hilbert transformation of the self energy 
  /// and chemical potential given as parameters.
  /// The density of states is specific for each of the derived classes
  ///
  /// @param G_omega the Greens function as a function of Matsubara Frequency omega 
  /// @param mu the chemical potential, h the magnetic field, beta the inverse temperature
  /// @return the result of the Hilbert transform: the bare Green's function G0 in Matsubara frequencies
  virtual matsubara_green_function_t operator()(const matsubara_green_function_t & G_omega, 
                                                matsubara_green_function_t &G0_omega, 
                                                const double mu, const double h, const double beta) const=0;
  virtual matsubara_green_function_t initial_G0(const alps::params& parms, const alps::params& input) const=0;
  
  template <class T>
  green_function<T> symmetrize(const green_function<T>& G, const bool symmetrization) const
  {
    green_function<T> G_new(G);
    if (symmetrization) {
      assert(G_new.nflavor()%2==0);
      for(spin_t flavor=0;flavor<G_new.nflavor(); flavor+=2){
        for(itime_index_t tau=0;tau<G_new.ntime();++tau){
          G_new(tau, flavor  )=0.5*(G_new(tau, flavor)+G_new(tau, flavor+1));
          G_new(tau, flavor+1)=G_new(tau, flavor);
        }
      }
    }
    return G_new;
  }

};


/// @brief performs a Hilbert transformation
///
/// The density of states is handled via class Bandstructure: 
///    semicircle DOS
///    user-defined via DOS histogram
///    tight-binding for square or hexagonal lattice
/// Arguments are expected to be in Matsubara Frequencies.
class GeneralFSHilbertTransformer : public FrequencySpaceHilbertTransformer {
public:
  
  GeneralFSHilbertTransformer(const alps::params& parms, const alps::params& input, bool /*ignored*/);
  GeneralFSHilbertTransformer(alps::params& parms, const alps::params& input);
  virtual ~GeneralFSHilbertTransformer() {}
  
  virtual matsubara_green_function_t operator()(const matsubara_green_function_t & G_omega, 
                                                matsubara_green_function_t &G0_omega, 
                                                const double mu, const double h, const double beta) const;
  virtual matsubara_green_function_t initial_G0(const alps::params& parms, const alps::params& input) const;

private:
  bool AFM;
  boost::shared_ptr<Bandstructure> bandstruct;
};


/// A Hilbert transformation for a semicircle density of states
/// Currently UNUSED: the FrequencySpaceHilbertTransformer is able to handle the semicircle DOS
/// NOTE: that the SemicircleFSHilbertTransformer::operator() uses the equation t^2*Delta=G for converged solution,
///       thus the result of the Hilbert transformation does not necesarilly equal to that of the 
///       FrequencySpaceHilbertTransformer::operator() in the not-yet-converged case
///       The effect on the convergency rate not fully explored.
class SemicircleFSHilbertTransformer : public FrequencySpaceHilbertTransformer {
public:
  explicit SemicircleFSHilbertTransformer(alps::params& parms)
    : bandstruct(parms) 
  {
    bandstruct.set_parms(parms);
  }

  /// It receives the dressed Green's function G(\omega) as input and returns tha bare GF G0(\omega)
  virtual matsubara_green_function_t operator()(const matsubara_green_function_t& G_omega, 
                                                matsubara_green_function_t &G0_omega_ignored, 
                                                const double mu, const double h, const double beta) const;
  
  virtual matsubara_green_function_t initial_G0(const alps::params& parms, const alps::params& input) const;
  
private:
  SemicircleBandstructure bandstruct;
};



#endif /*ALPS_DMFT_HILBERTTRANSFORMER_H*/
