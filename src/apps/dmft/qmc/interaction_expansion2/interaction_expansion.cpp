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

#include <boost/math/constants/constants.hpp>
#include "interaction_expansion.hpp"
#include <ctime>
#include "run_config.hpp"
#include "alps/ngs/make_deprecated_parameters.hpp"

//global variables

frequency_t c_or_cdagger::nm_;
bool c_or_cdagger::use_static_exp_;
unsigned int c_or_cdagger::ntau_;
double c_or_cdagger::beta_;
double *c_or_cdagger::omegan_;
std::complex<double> *c_or_cdagger::exp_iomegan_tau_;


InteractionExpansion::InteractionExpansion(const alps::run_configuration &run, int node)
: alps::mcbase(run.parameters,node),
max_order(run.parameters["MAX_ORDER"].as<unsigned int>()),
n_flavors(run.parameters["FLAVORS"].as<unsigned int>()),
n_site(run.parameters["SITES"].as<unsigned int>()),
n_matsubara(run.parameters["NMATSUBARA"].as<unsigned int>()),
n_matsubara_measurements(run.parameters["NMATSUBARA_MEASUREMENTS"].as<unsigned int>()),
n_tau(run.parameters["N"].as<unsigned int>()),
n_tau_inv(1./n_tau),
n_self(run.parameters["NSELF"].as<unsigned int>()),
mc_steps((boost::uint64_t)run.parameters["SWEEPS"]),
therm_steps(run.parameters["THERMALIZATION"].as<std::uint64_t>()),
beta((double)run.parameters["BETA"]),
temperature(1./beta),
onsite_U((double)run.parameters["U"]),
alpha((double)run.parameters["ALPHA"]),
U(alps::make_deprecated_parameters(run.parameters)),
recalc_period(run.parameters["RECALC_PERIOD"].as<unsigned int>()),
measurement_period(run.parameters["MEASUREMENT_PERIOD"].as<unsigned int>()),
almost_zero(run.parameters["ALMOSTZERO"].as<double>()),
green_matsubara(n_matsubara, n_site, n_flavors),
bare_green_matsubara(n_matsubara,n_site, n_flavors), 
bare_green_itime(n_tau+1, n_site, n_flavors),
green_itime(n_tau+1, n_site, n_flavors),
pert_hist(max_order)
{
  const auto &parms = run.parameters;
  random.engine().seed(run.execution["seed"].as<std::uint64_t>() + static_cast<std::uint64_t>(node));
  //initialize measurement method
  if (parms.value_or("HISTOGRAM_MEASUREMENT", false))
    measurement_method=selfenergy_measurement_itime_rs;
  else
    measurement_method=selfenergy_measurement_matsubara;
  for(unsigned int i=0;i<n_flavors;++i)
    g0.push_back(green_matrix(n_tau, 20));
  //other parameters
  weight=0;
  sign=1;
  step=0;
  measurement_time=0;
  update_time=0;
  read_ctint_bare_green(parms, run.input, bare_green_matsubara);
  if (run.input["atomic"].as<bool>()) {
    for (spin_t flavor = 0; flavor < n_flavors; ++flavor)
      for (itime_index_t t = 0; t <= n_tau; ++t)
        bare_green_itime(t, 0, 0, flavor) = -0.5;
  } else {
    FourierTransformer::generate_transformer(alps::make_deprecated_parameters(parms), fourier_ptr);
    fourier_ptr->backward_ft(bare_green_itime, bare_green_matsubara);
  }
  //initialize the simulation variables
  initialize_simulation(parms);
  if(node==0) {print(std::cout);}
  vertex_histograms=new simple_hist *[n_flavors*n_flavors];
  vertex_histogram_size=100;
  for(unsigned int i=0;i<n_flavors*n_flavors;++i){
    vertex_histograms[i]=new simple_hist(vertex_histogram_size);
  }
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
  if (is_thermalized())
    measure_observables();
}



double InteractionExpansion::fraction_completed() const{
  if (!is_thermalized())
    return 0.;
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
  nm_=p["NMATSUBARA"].as<unsigned int>();
  omegan_ = new double[nm_];
  for(unsigned int i=0;i<nm_;++i) {
    omegan_[i]=(2.*i+1.)*boost::math::constants::pi<double>()/beta_;
  }
  if(p.exists("TAU_DISCRETIZATION_FOR_EXP")) {
    ntau_=p["TAU_DISCRETIZATION_FOR_EXP"];
    use_static_exp_=true;
    exp_iomegan_tau_=new std::complex<double> [std::size_t(2)*nm_*ntau_];
    if(exp_iomegan_tau_==0){throw std::runtime_error("not enough memory for computing exp!"); }
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
