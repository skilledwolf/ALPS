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

#include <utility> // Boost.Math 1.76 includes this inside a namespace.
#include <boost/math/constants/constants.hpp>
#include "interaction_expansion.hpp"
#include "bare_green_input.hpp"
#include <ctime>
#include "xml.h"
#include "alps/ngs/make_deprecated_parameters.hpp"

//global variables

frequency_t c_or_cdagger::nm_;
bool c_or_cdagger::use_static_exp_;
unsigned int c_or_cdagger::ntau_;
double c_or_cdagger::beta_;
std::vector<double> c_or_cdagger::omegan_;
std::vector<std::complex<double>> c_or_cdagger::exp_iomegan_tau_;


InteractionExpansion::InteractionExpansion(const alps::params &parms, int node)
: alps::mcbase(parms,node),
max_order(parms["MAX_ORDER"] | 2048),
n_flavors(parms["FLAVORS"] | (parms["N_ORBITALS"] | 2)),
n_site(parms["SITES"] | 1),
n_matsubara((int)(parms["NMATSUBARA"]|parms["N_MATSUBARA"])),
n_matsubara_measurements(parms["NMATSUBARA_MEASUREMENTS"] | (int)n_matsubara),
n_tau((int)(parms["N"]|parms["N_TAU"])),
n_tau_inv(1./n_tau),
n_self(parms["NSELF"] | (int)(10*n_tau)),
mc_steps((boost::uint64_t)parms["SWEEPS"]),
therm_steps((unsigned int)parms["THERMALIZATION"]),        
max_time_in_seconds(parms["MAX_TIME"] | 86400),
beta((double)parms["BETA"]),                        
temperature(1./beta),
onsite_U((double)parms["U"]),                        
alpha((double)parms["ALPHA"]),
U(alps::make_deprecated_parameters(parms)),                         
recalc_period(parms["RECALC_PERIOD"] | 5000),
measurement_period(parms["MEASUREMENT_PERIOD"] | (parms["N_MEAS"] | 200)),
convergence_check_period(parms["CONVERGENCE_CHECK_PERIOD"] | (int)recalc_period),
almost_zero(parms["ALMOSTZERO"] | 1.e-16),
seed(parms["SEED"] | 0),
green_matsubara(n_matsubara, n_site, n_flavors),
bare_green_matsubara(n_matsubara,n_site, n_flavors), 
bare_green_itime(n_tau+1, n_site, n_flavors),
green_itime(n_tau+1, n_site, n_flavors),
pert_hist(max_order)
{
  //initialize measurement method
  if (parms["HISTOGRAM_MEASUREMENT"] | false)
    measurement_method=selfenergy_measurement_itime_rs;
  else
    measurement_method=selfenergy_measurement_matsubara;
  for(unsigned int i=0;i<n_flavors;++i)
    g0.push_back(green_matrix(n_tau, 20));
  //other parameters
  weight=0;
  sign=1;
  step=0;
  start_time=time(NULL);
  measurement_time=0;
  update_time=0;
  ctint_input::read_bare_green(parms,bare_green_matsubara);
  if(!(parms["ATOMIC"] | false)) {
    FourierTransformer::generate_transformer(alps::make_deprecated_parameters(parms), fourier_ptr);
    fourier_ptr->backward_ft(bare_green_itime, bare_green_matsubara);
  } else {
    for(spin_t flavor=0; flavor<n_flavors; ++flavor) 
      for(site_t site1=0; site1<n_site; ++site1) 
        for(site_t site2=0; site2<n_site; ++site2) 
          for(itime_index_t t=0; t<=n_tau; ++t)
            bare_green_itime(t, site1, site2, flavor)=-0.5;
  }
  //initialize the simulation variables
  initialize_simulation(parms);
  if(node==0) {print(std::cout);}
  c_or_cdagger::initialize_simulation(parms);
  
  if(n_site !=1) throw std::invalid_argument("you're trying to run this code for more than one site. Do you know what you're doing?!?");
}



void InteractionExpansion::update()
{
  for(std::size_t i=0;i<measurement_period;++i){
    step++;
    interaction_expansion_step();                
    if(vertices.size()<max_order)
      pert_hist[vertices.size()]++;
    if(step % recalc_period ==0)
      reset_perturbation_series();
  }
}
void InteractionExpansion::measure(){
  if (is_thermalized()) measure_observables();
}



double InteractionExpansion::fraction_completed() const{
  //check for error convergence
  if (!is_thermalized())
    return 0.;
  if(time(NULL)-start_time> max_time_in_seconds){
    std::cout<<"we ran out of time!"<<std::endl;
    return 1;
  }
  return ((step-therm_steps) / (double) mc_steps);
}



///do all the setup that has to be done before running the simulation.
void InteractionExpansion::initialize_simulation(const alps::params &parms)
{
  weight=0;
  sign=1;
  //set the right dimensions:
  for(spin_t flavor=0;flavor<n_flavors;++flavor)
    M.push_back(inverse_m_matrix());
  vertices.clear();
  pert_hist.clear();
  //initialize ALPS observables
  initialize_observables();
  green_matsubara=bare_green_matsubara;
  green_itime=bare_green_itime;
}



void c_or_cdagger::initialize_simulation(const alps::params &p)
{
  beta_=p["BETA"];
  // Adapted owning reinitialization from ALPS reconciliation
  // 90a36c0d4c64ee7308dcb0651264197ffdf36449, interaction_expansion2/interaction_expansion.cpp.
  // Operators compute all frequencies even when only a prefix is measured.
  // The former measured-prefix allocation was read past its end in compute_exp.
  nm_=p["NMATSUBARA"]|p["N_MATSUBARA"];
  omegan_.resize(nm_);
  exp_iomegan_tau_.clear();
  for(unsigned int i=0;i<nm_;++i) {
    omegan_[i]=(2.*i+1.)*boost::math::constants::pi<double>()/beta_;
  }
  if(p.defined("TAU_DISCRETIZATION_FOR_EXP")) {
    ntau_=p["TAU_DISCRETIZATION_FOR_EXP"];
    use_static_exp_=true;
    exp_iomegan_tau_.resize(std::size_t(2)*nm_*ntau_);
    std::cout<<"starting computation of exp values for measurement"<<std::endl;
    for(unsigned int i=0;i<ntau_;++i){
      double tau=i*beta_/(double)ntau_;
      for(unsigned int o=0;o<nm_;++o)
        exp_iomegan_tau_[std::size_t(2)*nm_*i + o] = std::complex<double>(cos(omegan_[o]*tau), sin(omegan_[o]*tau));
      for(unsigned int o=0;o<nm_;++o)
        exp_iomegan_tau_[std::size_t(2)*nm_*i + nm_ + o] = std::complex<double>(cos(omegan_[o]*tau), -sin(omegan_[o]*tau));
    }
    std::cout<<"done exp computation."<<std::endl;
  } else {
    use_static_exp_=false;
  }
}
