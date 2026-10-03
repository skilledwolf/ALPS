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
#include <alps/params.hpp>
#include <cmath>
#include <cassert>
#include <fstream>
#include <vector>

class U_matrix{
public:
  U_matrix(const alps::params &parms, const alps::params &input = {}) :
    ns_(parms.value_or<site_t>("SITES", 1)),
    nf_(parms.value_or<spin_t>("FLAVORS", 2)),
      n_nonzero_(0), mu_shift_(0)
  {
    if (ns_ != 1 || nf_ == 0) throw std::invalid_argument("Density interactions require SITES=1 and positive FLAVORS");
    val_.resize(std::size_t(nf_)*nf_, 0.); //default: non-interacting.
    if(input.exists("interaction_matrix")){
      const auto ufilename = input["interaction_matrix"].as<std::string>();
      std::ifstream u_file(ufilename.c_str());
      if (!u_file.is_open()) throw std::runtime_error("Cannot open input.interaction_matrix: " + ufilename);
      int i;
      int j;
      double U_ij;
      while (true) {
        u_file >> std::ws;
        if (u_file.eof()) break;
        if (!(u_file >> i >> j >> U_ij)) throw std::invalid_argument("Malformed input.interaction_matrix");
        if (i < 0 || j < 0 || static_cast<unsigned>(i) >= nf_ || static_cast<unsigned>(j) >= nf_ ||
            !std::isfinite(U_ij))
          throw std::invalid_argument("Invalid index or value in input.interaction_matrix");
        operator()(i,j)=U_ij;
      }
      if (u_file.bad()) throw std::runtime_error("Cannot read input.interaction_matrix");
    } else if (nf_==1) {
      //special case: only 1 orbital
      operator()(0,0)=parms["U"].as<double>();
    }else if (nf_==2) {
      //your ordinary two site problem
      assert(parms.exists("U"));
      double U=parms["U"].as<double>();
      operator()(0,0)=0; operator()(1,1)=0;
      operator()(0,1)=U; operator()(1,0)=U;
    } else {
      assert(parms.exists("U") && parms.exists("J"));
      double U=parms["U"].as<double>();
      double J=parms["J"].as<double>();
      double Uprime = parms.value_or("U'", U-2*J);
      assemble(U, Uprime, J);
    }
    for (std::size_t i=0; i<val_.size(); ++i)
      if (!std::isfinite(val_[i])) throw std::invalid_argument("Density interaction matrix must be finite");
      else if (val_[i]!=0)
        n_nonzero_++;
    for (unsigned i=0; i<nf_; ++i) 
      mu_shift_ += operator()(i,0);
    mu_shift_ /= 2;
    if (!std::isfinite(mu_shift_)) throw std::invalid_argument("Density interaction chemical-potential shift must be finite");
  }
  

  void assemble(const double U, const double Uprime, const double J){
    //this implements the U matrix for the special case of n_flavor/2 degenerate bands
    if (ns_ != 1 || nf_ % 2 != 0) throw std::invalid_argument("Hund interactions require SITES=1 and paired FLAVORS");
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
  std::vector<double> val_;
  site_t ns_;
  spin_t nf_;
  int n_nonzero_;
  double mu_shift_;
};

std::ostream &operator<<(std::ostream &os, const U_matrix &U);
//U_MATRIX_H
#endif 
