/****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2012 by Emanuel Gull <egull@umich.edu>,
 *                       Hartmut Hafermann <hafermann@cpht.polytechnique.fr>
 *
 *  based on an earlier version by Philipp Werner and Emanuel Gull
 *
 *
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
 *
 *****************************************************************************/
#ifndef HYB_SIM_MAIN
#define HYB_SIM_MAIN
#endif

#include <iomanip>
#include"hyb.hpp"

hybridization::hybridization(const alps::params &parms, int crank_)
: alps::mcbase(parms, crank_),
crank(crank_),
local_config(parms,crank),
hyb_config(parms)
{
  sanity_check(parms); //before doing anything, check whether the input parameters make sense
  show_info(parms,crank);
  
  //initializing general simulation constants
  nacc.assign(7,0.);
  nprop.assign(7,0.);
  sweep_count=0;
  output_period=parms.value_or("OUTPUT_PERIOD", 100000);
  //lasttime = boost::chrono::steady_clock::now();
  //delay = boost::chrono::seconds(parms.value_or("OUTPUT_PERIOD", 600));
  
  update_type.clear();
  update_type.push_back("change zero state   ");
  update_type.push_back("insert segment      ");
  update_type.push_back("remove segment      ");
  update_type.push_back("insert anti-segment ");
  update_type.push_back("remove anti-segment ");
  update_type.push_back("swap segment        ");
  update_type.push_back("global flip         ");
  
  sweeps=0;                                                                        //Sweeps currently done
  thermalization_sweeps = parms["THERMALIZATION"];                                 //Sweeps to be done for thermalization
  total_sweeps = parms["SWEEPS"];                                                  //Sweeps to be done in total
  n_orbitals = parms["N_ORBITALS"];                                                //number of orbitals
  sign = 1.;                                                                       //fermionic sign. plus or minus one.
  
  //initializing physics parameters
  beta = parms["BETA"];                                                            //inverse temperature
  
  //initializing updates parameters
  N_meas = parms["N_MEAS"];                                                        //number of updates per measurement
  N_hist_orders = parms.value_or("N_HISTOGRAM_ORDERS", 50);                                  //number of orders that are measured for the order histogram
  MEASURE_timeseries = parms.value_or("TIMESERIES", false);
  NUM_BINS = parms.value_or("NUM_BINS", 0);
  //std::cerr << "NUM_BINS = " << NUM_BINS << std::endl;
  //initializing measurement parameters
  spin_flip = parms.value_or("SPINFLIP", false);                                                //whether to perform local spin-flip updates
  global_flip = parms.value_or("GLOBALFLIP", false);                                                //whether to perform global spin-flip updates
  MEASURE_nnt = parms.value_or("MEASURE_nnt", false);                                           //measure density-density correlation function in imaginary time
  MEASURE_nnw = parms.value_or("MEASURE_nnw", false);                                           //measure density-density correlation function in frequency
  MEASURE_nn = parms.value_or("MEASURE_nn", false);                                              //measure density-density correlation function at equal times
  MEASURE_g2w = parms.value_or("MEASURE_g2w", false);                                            //measure two-particle Green function
  MEASURE_h2w = parms.value_or("MEASURE_h2w", false);                                            //measure higher-order two-particle correlator
  MEASURE_time = parms.value_or("MEASURE_time", true);                                          //measure in imaginary time (ON by default)
  MEASURE_freq = parms.value_or("MEASURE_freq", false);                                          //measure in frequency space
  MEASURE_legendre = parms.value_or("MEASURE_legendre", false);                                  //measure in legendre polynomials
  MEASURE_sector_statistics = parms.value_or("MEASURE_sector_statistics", false);                //measure sector statistics
  N_w = parms.value_or("N_MATSUBARA", 0);                                                    //number of Matsubara frequencies for gw
  N_l = parms.value_or("N_LEGENDRE", 0);                                                     //number of Legendre polynomial coefficients
  N_t = parms["N_TAU"];                                                            //number of tau slices for gt
  N_nn = parms.value_or("N_nn", 0);                                                          //number of tau-points on which density density correlator is measured
  N_w2 = parms.value_or("N_w2", 0);                                                          //number of Matsubara frequency points for two-particle measurements
  N_W = parms.value_or("N_W", 0);                                                            //number of bosonic Matsubara frequency points for two-particle measurements
  N_w_aux = (N_w2+N_W>1 ? N_w2+N_W-1 : 0);                                         //number of Matsubara frequency points for the measurment of M(w1,w2)
  
  //create measurement objects
  create_measurements();
  
  if(crank==0){
    std::cout<<"Hybridization Expansion Simulation CT-HYB"<<std::endl;
    std::cout<<"Part of the ALPS DMFT Project"<<std::endl;
    std::cout<<"Refer to the documentation for more information."<<std::endl;
  }
  
  start_time=clock();
  end_time=start_time+ CLOCKS_PER_SEC*((long)parms["MAX_TIME"]);

  
  //std::cout<<"process " << crank << " starting simulation"<<std::endl;
  csize=1;
 //we don't have a nice way of getting the MPI size from ALPS, because we don't know about the communicator at this point.
 //here is a safe way of getting the pool size into csize.
#ifdef ALPS_HAVE_MPI
  int mpi_init;
  MPI_Initialized(&mpi_init);
  if(mpi_init){
     MPI_Comm_size(MPI_COMM_WORLD, &csize);
  }
#endif
  std::cout<<"process " << crank << " of total: "<<csize<<" starting simulation"<<std::endl;
}

void hybridization::sanity_check(const alps::params &parms){
  //check whether the input parameters make sense before computing
  //NOTE: these checks are likely not to be complete, passing all checks does not guarantee all parameters to be meaningful!
  
  //first check that all mandatory parameters are defined
  if(!parms.exists("N_TAU")) throw std::invalid_argument("please specify the parameter N_TAU");
  if(!parms.exists("BETA")) throw std::invalid_argument("please specify parameter BETA for inverse temperature");
  if(!parms.exists("N_MEAS")) throw std::invalid_argument("please specify parameter N_MEAS for measurement interval");
  if(!parms.exists("THERMALIZATION") ||
     !parms.exists("SWEEPS") ||
     !parms.exists("N_ORBITALS") ) throw std::invalid_argument("please specify parameters THERMALIZATION, SWEEPS, and N_ORBITALS");
  
  //check paramater that are conditionally required
  if(parms.value_or("MEASURE_freq", false) && !parms.exists("N_MATSUBARA")) throw std::invalid_argument("please specify parameter N_MATSUBARA for # of Matsubara frequencies to be measured");
  
  if(parms.value_or("MEASURE_legendre", false) && !parms.exists("N_LEGENDRE")) throw std::invalid_argument("please specify parameter N_LEGENDRE for # of Legendre coefficients to be measured");
  if(parms.value_or("MEASURE_legendre", false) && !parms.exists("N_MATSUBARA")) throw std::invalid_argument("please specify parameter N_MATSUBARA for # of Matsubara frequencies");
  if(parms.value_or("MEASURE_nnt", false) && !parms.exists("N_nn")) throw std::invalid_argument("please specify the parameter N_nn for # of imaginary time points for the density-density correlator");
  if(parms.value_or("MEASURE_nnw", false) && !parms.exists("N_W")) throw std::invalid_argument("please specify the parameter N_W for # of bosonic frequencies for the density-density correlator");
  if(parms.value_or("MEASURE_g2w", false) || parms.value_or("MEASURE_h2w", false) ){
    if(!parms.exists("N_w2") ) throw std::invalid_argument("please specify the parameter N_w2 for # of fermionic Matsubara frequencies for two-particle functions");
    if(!parms.exists("N_W") ) throw std::invalid_argument("please specify the parameter N_W for # of bosonic Matsubara frequencies for two-particle functions");
    if((int)parms["N_w2"]%2!=0) throw std::invalid_argument("parameter N_w2 must be even");
  }
  if(parms.value_or("COMPUTE_VERTEX", false)){
    if( !(parms.value_or("MEASURE_freq", false)) ) throw std::invalid_argument("frequency measurement is required for computing the vertex, please set MEASURE_freq=1");
    
    if(! (parms.value_or("MEASURE_g2w", false) || parms.value_or("MEASURE_h2w", false) ) ) throw std::invalid_argument("at least one two-particle quantity is required for computing the vertex, set MEASURE_g2w=1 or MEASURE_h2w=1");
    if((int) parms["N_MATSUBARA"] < ((int)parms["N_w2"]/2 + (int)parms["N_W"] - 1) ) throw std::invalid_argument("for computing the vertex, N_MATSUBARA must be at least N_w2/2+N_W-1");
  }
  VERBOSE = (parms.value_or("VERBOSE", false));
  
  return;
}


void hybridization::show_info(const alps::params &parms, int crank){
  if(!(parms.value_or("VERBOSE", false))) return;

  //provide info on what is measured and how long the simulation will run
  if(!crank){
    if(parms.value_or("MEASURE_time", true)) std::cout << "measuring gt" << std::endl;
    if(parms.value_or("MEASURE_freq", false)) std::cout << "measuring gw" << std::endl << "measuring fw" << std::endl;
    if(parms.value_or("MEASURE_legendre", false)) std::cout << "measuring gl" << std::endl << "measuring fl" << std::endl;
    if(parms.value_or("MEASURE_g2w", false)) std::cout << "measuring g2w" << std::endl;
    if(parms.value_or("MEASURE_h2w", false)) std::cout << "measuring h2w" << std::endl;
    if(parms.value_or("MEASURE_nn", false)) std::cout << "measuring nn" << std::endl;
    if(parms.value_or("MEASURE_nnt", false)) std::cout << "measuring nnt" << std::endl;
    if(parms.value_or("MEASURE_nnw", false)) std::cout << "measuring nnw" << std::endl;
    if(parms.value_or("MEASURE_sector_statistics", false)) std::cout << "measuring sector statistics" << std::endl;
    if(parms.value_or("COMPUTE_VERTEX", false)) std::cout << "vertex will be computed" << std::endl;
    if(parms.exists("RET_INT_K")) std::cout << "using retarded interaction" << std::endl;
    if(parms.exists("U_MATRIX")) std::cout << "reading U matrix from file " << parms["U_MATRIX"] << std::endl;
    if(parms.exists("MU_VECTOR")) std::cout << "reading MU vector from file " << parms["MU_VECTOR"] << std::endl;
    std::cout << "Simulation scheduled to run " << parms["MAX_TIME"] << " seconds" << std::endl << std::endl;
  }
  return;
}

//this is a debug function that recomputes the full weight of the local and the hybridization configuration
double hybridization::full_weight() const{
  return local_config.full_weight()*hyb_config.full_weight();
}

std::ostream &operator<<(std::ostream &os, const hybridization &hyb){
  os<<cred<<"-----------------------------------------------------------------------------------"<<cblack<<std::endl;
  os<<hyb.local_config<<std::endl;
  os<<hyb.hyb_config<<std::endl;
  os<<cred<<"-----------------------------------------------------------------------------------"<<cblack<<std::endl;
  return os;
}
std::ostream &operator<<(std::ostream &os, const segment &s){
  os<<"( "<<s.t_start_<<" , "<<s.t_end_<<" ) ";
  return os;
}
double hybridization::fraction_completed()const{
  if(!is_thermalized()) return 0.;
  double work_fraction= (sweeps-thermalization_sweeps)/(double)total_sweeps;
  double time_fraction= (clock()-start_time)/(double)(end_time-start_time);
  //return max of sweeps done and time used. Divide time used by the number of processes in pool (all work done will be added up)
  return std::max(work_fraction, time_fraction/csize);
  //return work_fraction;
}
