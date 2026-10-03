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

#include"hybfun.hpp"
#include "input.hpp"
//construct a hybridization function. ntime: number of time slices.
//noffidag_orbitals: number of offdiagonal orbitals. ndiag_orbitals: number
//of diagonal orbitals. 
hybfun::hybfun(const alps::params &p, const alps::params &input):
green_function<double>(p["N_TAU"].as<int>()+1, 1, p["N_ORBITALS"])
{
  if(!p.exists("N_TAU") || (int)(p["N_TAU"])==0) throw std::invalid_argument("define parameter N_TAU, the number of hybridization time slices!");
  beta_=p["BETA"];

  //read in Green's function from a file
  read_hybridization_function(p,input);
  hybridization_function_sanity_check();

  extern int global_mpi_rank;
  if(global_mpi_rank==0){
    std::cout<<*this<<std::endl;
  }
}

void hybfun::hybridization_function_sanity_check(void){
for(std::size_t i=0; i<ntime();++i)
  for(std::size_t j=0; j<nflavor();++j)
    if(operator()(i,j)>0.) {
      std::cerr << "ERROR: Delta(t="<<i<<"; f="<<j<<") = " << operator()(i,j) << "  is positive." << std::endl;
      std::cerr << "Note: small positive values might be due to noise, in that case try to increase execution.time_limit or parameters.SWEEPS." << std::endl << std::flush;
      throw std::invalid_argument("Problem with hybridization function: Delta(\\tau) > 0. Delta should always be negative!");
    }
}

//this routine reads in the hybridization function, either from a text file or from an hdf5 file (for easy passing of binary data).
//In case of text files the file format is index - hyb_1 - hyb2 - hyb3 - ... in columns that go from time=0 to time=beta. Note that
//the hybridization function is in imaginary time and always positive between zero and \beta.
void hybfun::read_hybridization_function(const alps::params &p, const alps::params &input){
  const auto values=cthyb_input::series(p,input,false);
  for(std::size_t i=0;i<ntime();++i)
    for(std::size_t j=0;j<nflavor();++j) operator()(i,j)=values[i*nflavor()+j];
}

std::ostream &operator<<(std::ostream &os, const hybfun &hyb){
  os<<"the hybridization function is: "<<std::endl;
  for(std::size_t i=0;i<std::min<std::size_t>(10,hyb.ntime());++i){ std::cout<<i<<" "; for(std::size_t j=0;j<hyb.nflavor();++j){ std::cout<<hyb(i,j)<<" ";} std::cout<<std::endl; }
  os<<"... *** etc *** ...\n";
  os<<hyb.ntime()-1<<" "; for(std::size_t j=0;j<hyb.nflavor();++j){ std::cout<<hyb(hyb.ntime()-1,j)<<" ";} std::cout<<std::endl;
  return os;
}

//linear interpolation of the hybridization function. 
double hybfun::interpolate(double time, int orbital) const{

  ///TODO: do a spline here, write a better/faster interpolation routine.
  double sign=1;
  if (time<0) {
    time += beta_;
    sign=-1;
  }

  //this is the overall flip of Delta (in comparison to F)
  sign*=-1;
  //the code takes Delta as input, but internally works with F(tau)=-Delta(beta-tau)
  time=beta_-time;

  double n = time/beta_*(ntime()-1);
  int n_lower = (int)n; // interpolate linearly between n_lower and n_lower+1

  return sign*(operator()(n_lower,orbital) + (n-n_lower)*(operator()(n_lower+1,orbital)-operator()(n_lower,orbital)));
}
