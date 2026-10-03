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

#include"hybretintfun.hpp"
#include "input.hpp"
//construct a retarded interaction function. ntime: number of time slices.
//noffidag_orbitals: number of offdiagonal orbitals. ndiag_orbitals: number
//of diagonal orbitals. 
ret_int_fun::ret_int_fun(const alps::params &p, const alps::params &input):
green_function<double>(p["N_TAU"].as<int>()+1, 1, 2)//2 "flavors": K and K'
{
  bool use_retarded_interaction=input.exists("retarded_interaction");
  if(!use_retarded_interaction) return;

  if(!p.exists("N_TAU") || (int)(p["N_TAU"])==0) throw std::invalid_argument("define parameter N_TAU, the number of retarded interaction time slices!");
  beta_=p["BETA"];
  
  //read in Green's function from a file
  read_interaction_K_function(p,input);
  interaction_K_function_sanity_check();

  extern int global_mpi_rank;
  if(global_mpi_rank==0){
    std::cout<<*this<<std::endl;
  }
}

void ret_int_fun::interaction_K_function_sanity_check(void){
for(std::size_t i=0; i<ntime();++i)
  for(std::size_t j=0; j<nflavor();++j)
    if(operator()(i,j) < 0. && j==0) throw std::invalid_argument("Problem with retarded interaction function: RET_INT_K(\\tau) < 0. K should always be positive!");
if(operator()(0,0)!=0.) throw std::invalid_argument("Problem with retarded interaction function: RET_INT_K(\\tau=0) must be zero.");
}

//this routine reads in the retarded interaction function, either from a text file or from an hdf5 file (for easy passing of binary data).
//In  case of text files the file format is index - hyb_1 - hyb2 - hyb3 - ... in columns that go from time=0 to time=beta. Note that
//the retarded interaction function is in imaginary time and always positive both for negative and positive times. It is also symmetric.
void ret_int_fun::read_interaction_K_function(const alps::params &p, const alps::params &input){
  const auto values=cthyb_input::series(p,input,true);
  for(std::size_t i=0;i<ntime();++i)
    for(std::size_t j=0;j<nflavor();++j) operator()(i,j)=values[i*nflavor()+j];
}

std::ostream &operator<<(std::ostream &os, const ret_int_fun &K){
  os<<"the retarded interaction function and derivative are: "<<std::endl;
  for(std::size_t i=0;i<std::min<std::size_t>(10,K.ntime());++i){ std::cout<<i<<" "; for(std::size_t j=0;j<K.nflavor();++j){ std::cout<<K(i,j)<<" ";} std::cout<<std::endl; }
  os<<"... *** etc *** ...\n";
  os<<K.ntime()-1<<" "; for(std::size_t j=0;j<K.nflavor();++j){ std::cout<<K(K.ntime()-1,j)<<" ";} std::cout<<std::endl;
  return os;
}
//linear interpolation of the retarded interaction function 
double ret_int_fun::interpolate(double time) const{
  time=std::abs(time);

  double n = time/beta_*(ntime()-1);
  int n_lower = (int)n; // interpolate linearly between n_lower and n_lower+1

  return operator()(n_lower,0) + (n-n_lower)*(operator()(n_lower+1,0)-operator()(n_lower,0));
}

//linear interpolation of the retarded interaction function.
double ret_int_fun::interpolate_deriv(double time) const{
  //if(time<-beta_ || time > beta_) std::cout<< "!!" << std::endl;
  if(time<0.) return -1.*interpolate_deriv(-time);

  double n = time/beta_*(ntime()-1);
  int n_lower = (int)n; // interpolate linearly between n_lower and n_lower+1

  return operator()(n_lower,1) + (n-n_lower)*(operator()(n_lower+1,1)-operator()(n_lower,1));
}
