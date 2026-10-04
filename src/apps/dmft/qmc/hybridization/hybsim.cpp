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

hybridization::hybridization(const alps::run_configuration &run, int crank_)
: verbose(run.parameters["VERBOSE"].as<bool>()),
crank(crank_),
local_config(run.parameters,run.input,crank),
hyb_config(run.parameters,run.input)
{
  const auto &parms=run.parameters;
  random.engine().seed(run.execution["seed"].as<std::uint32_t>()+static_cast<std::uint32_t>(crank));
  show_info(run,crank);
  
  //initializing general simulation constants
  nacc.assign(7,0.);
  nprop.assign(7,0.);
  sweep_count=0;
  output_period=run.execution["progress_period"].as<int>();
  
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
  create_measurements(run.execution["bins"].as<std::size_t>());
  
  if(crank==0){
    std::cout<<"Hybridization Expansion Simulation CT-HYB"<<std::endl;
    std::cout<<"Part of the ALPS DMFT Project"<<std::endl;
    std::cout<<"Refer to the documentation for more information."<<std::endl;
  }
  
  std::cout<<"process " << crank << " starting simulation"<<std::endl;
}

void hybridization::run(std::function<bool()> const& stop_callback) {
  while (!stop_callback() && fraction_completed() < 1.) {
    update();
    measure();
  }
}

void hybridization::show_info(const alps::run_configuration &run, int crank){
  const auto &parms=run.parameters;
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
    if(run.input.exists("retarded_interaction")) std::cout << "using retarded interaction" << std::endl;
    if(run.input.exists("interaction_matrix")) std::cout << "reading U matrix from file " << run.input["interaction_matrix"] << std::endl;
    if(run.input.exists("chemical_potential")) std::cout << "reading MU vector from file " << run.input["chemical_potential"] << std::endl;
    std::cout << "Simulation scheduled to run " << run.execution["time_limit"] << " seconds" << std::endl << std::endl;
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
  return (sweeps-thermalization_sweeps)/(double)total_sweeps;
}
