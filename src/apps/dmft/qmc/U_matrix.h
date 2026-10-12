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


//Data structure for repulsion for density density bands.
#ifndef U_MATRIX_H
#define U_MATRIX_H
#include "types.h"
#include "alps/parameter.h"
#include <cmath>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <vector>

class U_matrix{
public:
  // Owned storage and checked input adapt ALPS reconciliation
  // 90a36c0d4c64ee7308dcb0651264197ffdf36449, this header. Retain the existing parameter/file interface;
  // vector storage makes copied matrices independent and exception-safe.
  U_matrix(const alps::Parameters &parms) :
    ns_(dimension(parms,"SITES",1)), nf_(dimension(parms,"FLAVORS",2)),
    n_nonzero_(0), mu_shift_(0) {
    if (nf_>static_cast<unsigned>(std::numeric_limits<int>::max())/nf_)
      throw std::invalid_argument("Density interaction dimensions exceed the retained count range");
    const auto cells=std::size_t(nf_)*nf_;
    if (cells>static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
        cells>std::numeric_limits<std::size_t>::max()/sizeof(double))
      throw std::invalid_argument("Density interaction matrix dimensions exceed the retained count/storage range");
    val_.assign(cells,0.);
    if (parms.defined("U_MATRIX")) {
      const std::string filename(parms["U_MATRIX"]);
      std::ifstream file(filename);
      if (!file) throw std::runtime_error("Cannot open density interaction matrix: "+filename);
      while (true) {
        file >> std::ws;
        if (file.eof()) break;
        int first,second; double interaction;
        if (!(file>>first>>second>>interaction)) throw std::invalid_argument("Malformed density interaction matrix triple");
        if (first<0 || second<0 || static_cast<unsigned>(first)>=nf_ || static_cast<unsigned>(second)>=nf_ || !std::isfinite(interaction))
          throw std::invalid_argument("Density interaction index or value is invalid");
        operator()(first,second)=interaction;
      }
      if (file.bad()) throw std::runtime_error("Cannot read density interaction matrix");
    } else if (nf_==1) operator()(0,0)=static_cast<double>(parms["U"]); // Retained half-filling representation.
    else if (nf_==2) operator()(0,1)=operator()(1,0)=static_cast<double>(parms["U"]);
    else {
      const double interaction=parms["U"],hund=parms.value_or_default("J",0.);
      const double interband=parms.value_or_default("U'",interaction-2.*hund);
      assemble(interaction,interband,hund);
    }
    for (spin_t first=0; first<nf_; ++first) {
      if (nf_>1 && operator()(first,first)!=0.) throw std::invalid_argument("Density interaction diagonal must vanish");
      for (spin_t second=0; second<nf_; ++second) {
        const double interaction=operator()(first,second);
        if (!std::isfinite(interaction) || (nf_>1 && interaction!=operator()(second,first)))
          throw std::invalid_argument("Density interaction matrix must be finite and symmetric");
        n_nonzero_+=interaction!=0.;
      }
      mu_shift_+=operator()(first,0);
    }
    mu_shift_/=2.;
    if (!std::isfinite(mu_shift_)) throw std::invalid_argument("Density interaction shift is not representable");
  }

  

  void assemble(const double U, const double Uprime, const double J){
    //this implements the U matrix for the special case of n_flavor/2 degenerate bands
    if (ns_!=1 || nf_%2!=0) throw std::invalid_argument("Generated Hund interactions require one site and paired flavors");
    for(spin_t i=0;i<nf_;i+=2){
      operator()(i  , i  ) = 0; //Pauli
      operator()(i+1, i+1) = 0; //Pauli
      operator()(i  , i+1) = U; //Hubbard repulsion same band
      operator()(i+1, i  ) = U; //Hubbard repulsion same band
      for(spin_t j=0; j<nf_; j+=2){
        if(j==i) 
          continue;
        operator()(i  ,j  ) = Uprime-J; //Hubbard repulsion interband same spin 
        operator()(i+1,j+1) = Uprime-J; //Hubbard repulsion interband same spin
        operator()(i  ,j+1) = Uprime; //Hubbard repulsion interband opposite spin (this used to be '+J', the rest of the world uses '-J' -> changed to be consistent).
        operator()(i+1,j  ) = Uprime; //Hubbard repulsion interband opposite spin
      }
    }
  } 
  
  const std::vector<double>& data() const { return val_; }
  
  double &operator()(spin_t flavor_i, spin_t flavor_j){
    return val_[flavor_i*nf_+flavor_j];
    }
  
  const double &operator() (spin_t flavor_i, spin_t flavor_j)const {
    return val_[flavor_i*nf_+flavor_j];
  }
  
  spin_t nf()const {return nf_;}
  spin_t ns()const {return ns_;}
  double mu_shift() const { return mu_shift_; }

  inline int n_nonzero() const{return n_nonzero_;}

private:
  static unsigned dimension(const alps::Parameters& parms,const char* key,int fallback) {
    const long long value=parms.defined(key) ? static_cast<long long>(parms[key]) : fallback;
    if (value<1 || static_cast<unsigned long long>(value)>std::numeric_limits<unsigned>::max())
      throw std::invalid_argument("Density interaction dimensions must be positive integers in range");
    return static_cast<unsigned>(value);
  }
  std::vector<double> val_;
  site_t ns_;
  spin_t nf_;
  int n_nonzero_;
  double mu_shift_;
};

std::ostream &operator<<(std::ostream &os, const U_matrix &U);
//U_MATRIX_H
#endif 
