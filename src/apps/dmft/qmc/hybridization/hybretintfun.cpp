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
#include <algorithm>
#include <cmath>
#include <stdexcept>
//construct a retarded interaction function. ntime: number of time slices.
//noffidag_orbitals: number of offdiagonal orbitals. ndiag_orbitals: number
//of diagonal orbitals. 
ret_int_fun::ret_int_fun(const alps::params &p):
green_function<double>(p["N_TAU"]+1, 1, 2)//2 "flavors": K and K'
{
  bool use_retarded_interaction=p.defined("RET_INT_K");
  if(!use_retarded_interaction) return;

  if(!p.defined("N_TAU") || (int)(p["N_TAU"])==0) throw std::invalid_argument("define parameter N_TAU, the number of retarded interaction time slices!");
  beta_=p["BETA"];
  if (!std::isfinite(beta_) || beta_<=0.)
    throw std::invalid_argument("Retarded table BETA must be finite and positive");
  
  //read in Green's function from a file
  read_interaction_K_function(p);
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
void ret_int_fun::read_interaction_K_function(const alps::params &p){
  if(!p.defined("RET_INT_K")) throw(std::invalid_argument(std::string("Parameter RET_INT_K missing, filename for retarded interaction function not specified.")));
  std::string fname=p["RET_INT_K"].cast<std::string>();
  if(p.defined("K_IN_HDF5") && p["K_IN_HDF5"].cast<bool>()){//attempt to read from h5 archive
    alps::hdf5::archive ar(fname, alps::hdf5::archive::READ);
    std::vector<double> tmp(ntime());
    ar>>alps::make_pvp("/Ret_int_K",tmp);
      for(std::size_t i=0; i<ntime(); i++)
        operator()(i,0)=tmp[i];
    tmp.clear();
    ar>>alps::make_pvp("/Ret_int_Kp",tmp);
      for(std::size_t i=0; i<ntime(); i++)
        operator()(i,1)=tmp[i];
    tmp.clear();
  }
  else{//read from text file
    std::ifstream infile(fname.c_str());
    if(!infile.good()){
      throw(std::invalid_argument(std::string("could not open retarded interaction file (text format) ") + p["RET_INT_K"]));
    }
    for (std::size_t i=0; i<ntime(); i++) {
      if(!infile.good()){
        throw(std::invalid_argument(std::string("could not read retarded interaction file (text format) ") + p["RET_INT_K"] + ". probably wrong number of lines"));
      }
      double dummy;
      infile >> dummy;
      //std::cout<<i<<" "<<ntime()<<" "<<"dummy: "<<dummy<<std::endl;
      for (std::size_t j=0; j<nflavor(); j++){
        if(!infile.good()){
          throw(std::invalid_argument(std::string("could not read retarded interaction file (text format) ") + p["RET_INT_K"] + ". probably wrong number of columns"));
        }
        double delta;
        infile >> delta;
        operator()(i,j)=delta;
      }
    }
  }
}

std::ostream &operator<<(std::ostream &os, const ret_int_fun &K){
  // Short-mesh bound adapted from ALPS reconciliation
  // 90a36c0d4c64ee7308dcb0651264197ffdf36449, hybridization/hybretintfun.cpp.
  os<<"the retarded interaction function and derivative are: "<<std::endl;
  for(std::size_t i=0;i<std::min<std::size_t>(10,K.ntime());++i){ os<<i<<" "; for(std::size_t j=0;j<K.nflavor();++j){ os<<K(i,j)<<" ";} os<<std::endl; }
  os<<"... *** etc *** ...\n";
  os<<K.ntime()-1<<" "; for(std::size_t j=0;j<K.nflavor();++j){ os<<K(K.ntime()-1,j)<<" ";} os<<std::endl;
  return os;
}
//linear interpolation of the retarded interaction function 
double ret_int_fun::interpolate(double time) const{
  if (!std::isfinite(time) || std::abs(time) > beta_)
    throw std::invalid_argument("Retarded interpolation time must lie in [-BETA,BETA]");
  time=std::abs(time);

  double n = time/beta_*(ntime()-1);
  if (n >= ntime()-1) return operator()(ntime()-1,0);
  int n_lower = (int)n; // interpolate linearly between n_lower and n_lower+1

  return operator()(n_lower,0) + (n-n_lower)*(operator()(n_lower+1,0)-operator()(n_lower,0));
}

//linear interpolation of the retarded interaction function.
double ret_int_fun::interpolate_deriv(double time) const{
  if (!std::isfinite(time) || std::abs(time) > beta_)
    throw std::invalid_argument("Retarded derivative time must lie in [-BETA,BETA]");
  //if(time<-beta_ || time > beta_) std::cout<< "!!" << std::endl;
  if(time<0.) return -1.*interpolate_deriv(-time);

  double n = time/beta_*(ntime()-1);
  if (n >= ntime()-1) return operator()(ntime()-1,1);
  int n_lower = (int)n; // interpolate linearly between n_lower and n_lower+1

  return operator()(n_lower,1) + (n-n_lower)*(operator()(n_lower+1,1)-operator()(n_lower,1));
}
