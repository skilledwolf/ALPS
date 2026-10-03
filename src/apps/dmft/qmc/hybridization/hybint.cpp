/****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2012 by Emanuel Gull <gull@pks.mpg.de>,
 *                   
 *  based on an earlier version by Philipp Werner and Emanuel Gull
 *
 *
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
 *
 *****************************************************************************/

#include"hybint.hpp"
#include "input.hpp"


interaction_matrix::interaction_matrix(const alps::params &p, const alps::params &input){
  extern int global_mpi_rank;
  n_orbitals_=p["N_ORBITALS"];
  val_.resize(n_orbitals_*n_orbitals_,0.);
  if(input.exists("interaction_matrix")){
    val_=cthyb_input::static_values(input,"interaction_matrix","interaction_format","/Umatrix",n_orbitals_*n_orbitals_);
  }else{
    if(!p.exists("U")) throw std::invalid_argument("Specify parameters.U (optionally U' and J) or input.interaction_matrix");
    double U=(double)(p["U"]);
    double J=(double)(p.value_or("J", 0.));
    double Uprime = (p.value_or("U'", (U-2*J)));
    assemble(U, Uprime, J);
  }
}

void interaction_matrix::apply_shift(const double shift){
  for(int i=0; i<n_orbitals(); ++i)
    for(int j=0; j<n_orbitals(); ++j)
      if(i!=j) operator()(i,j)+=shift; //apply shift
}

void interaction_matrix::assemble(const double U, const double Uprime, const double J){
  if(Uprime==U && J==0){
     for(int i=0;i<n_orbitals_;++i){
       for(int j=0;j<n_orbitals_;++j){
         operator()(i,j)=(i==j)?0:U;
       }
     }
   }else{
  if(n_orbitals_%2!=0){
    std::cerr<<"n_orbitals is: "<<n_orbitals_<<std::endl;
    throw std::logic_error("extend assemble or write interaction matrix to file for odd # orbitals");
  }
  for(int i=0;i<n_orbitals_;i+=2){
    operator()(i  , i  ) = 0; //Pauli
    operator()(i+1, i+1) = 0; //Pauli
    operator()(i  , i+1) = U; //Hubbard repulsion same band
    operator()(i+1, i  ) = U; //Hubbard repulsion same band
    for(int j=0; j<n_orbitals_; j+=2){
      if(j==i) 
        continue;
      operator()(i  ,j  ) = Uprime-J; //Hubbard repulsion interband same spin 
      operator()(i+1,j+1) = Uprime-J; //Hubbard repulsion interband same spin
      operator()(i  ,j+1) = Uprime; //Hubbard repulsion interband opposite spin (this used to be '+J', the rest of the world uses '-J' -> changed to be consistent).
      operator()(i+1,j  ) = Uprime; //Hubbard repulsion interband opposite spin
    }
  }
}
} 

std::ostream &operator<<(std::ostream &os, const interaction_matrix &U){
  os<<"(effective) U matrix with "<<U.n_orbitals()<<" orbitals: "<<std::endl;
  for(int i=0;i<U.n_orbitals();++i){
    for(int j=0;j<U.n_orbitals();++j){
      os<<U(i,j)<<" ";
    }
    os<<std::endl;
  }
  return os;
}
std::ostream &operator<<(std::ostream &os, const chemical_potential &mu){
  os<<"(effective) chemical potential with "<<mu.n_orbitals()<<" orbitals: "<<std::endl;
  for(std::size_t i=0;i<mu.n_orbitals();++i){
    os<<mu[i]<<" ";
  }
  os<<std::endl;
  return os;
}

