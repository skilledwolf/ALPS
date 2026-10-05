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

#include "interaction_expansion.hpp"
#include "run_config.hpp"
#include <complex>

void evaluate_selfenergy_measurement_matsubara(const InteractionExpansion::results_type &results,
                                                                        matsubara_green_function_t &green_matsubara_measured,
                                                                        const matsubara_green_function_t &bare_green_matsubara, 
                                                                        const double &beta, std::size_t n_site, 
                                                                        std::size_t n_flavors, std::size_t n_matsubara);
void evaluate_selfenergy_measurement_itime_rs(const InteractionExpansion::results_type &results,
                                                                       itime_green_function_t &green_result,
                                                                       const itime_green_function_t &green0, 
                                                                       const double &beta, const int n_site, 
                                                                       const int n_flavors, const int n_tau, const int n_self);



void compute_greens_functions(const InteractionExpansion::results_type &results, const alps::params& parms, const alps::params &input, alps::hdf5::archive& archive)
{
  std::cout<<"getting result!"<<std::endl;
  unsigned int n_matsubara = parms["NMATSUBARA"].as<unsigned int>();
  unsigned int n_matsubara_measurements=parms["NMATSUBARA_MEASUREMENTS"].as<unsigned int>();
  unsigned int n_tau=parms["N"].as<unsigned int>();
  unsigned int n_self=parms["NSELF"].as<unsigned int>();
  spin_t n_flavors(parms["FLAVORS"].as<unsigned int>());
  unsigned int n_site(parms["SITES"].as<unsigned int>());
  double beta(parms["BETA"]);
  itime_green_function_t green_itime_measured(n_tau+1, n_site, n_flavors);
  matsubara_green_function_t green_matsubara_measured(n_matsubara, n_site, n_flavors);
  boost::shared_ptr<FourierTransformer> fourier_ptr;
  boost::shared_ptr<FourierTransformer> fourier_ptr_g0;
  FourierTransformer::generate_transformer(parms, fourier_ptr_g0);
  //find whether our data is in imaginary time or frequency:
  bool measure_in_matsubara=true;
  if(parms.value_or("HISTOGRAM_MEASUREMENT", false))
    measure_in_matsubara=false;
  auto mean_order = results.at("PertOrder").mean();
  
  std::cout<<"average matrix size was: "<<std::endl;
  for(unsigned int i=0;i<n_flavors;++i){
    std::cout<<mean_order[i]<<"\t";
  }
  std::cout<<std::endl;
  std::cout<<"average sign was: "<<results.at("Sign").mean()(0)<<" error: "<<results.at("Sign").stderror()(0)<<std::endl;
  //single particle Green function measurements
  matsubara_green_function_t bare_green_matsubara(n_matsubara, n_site, n_flavors);
  const auto density_mean = results.at("densities").mean();
  const std::vector<double> densities(density_mean.data(), density_mean.data()+density_mean.size());
  read_ctint_bare_green(parms, input, bare_green_matsubara);
  if(measure_in_matsubara) {
    evaluate_selfenergy_measurement_matsubara(results, green_matsubara_measured, 
                                              bare_green_matsubara,
                                              beta, n_site, n_flavors, n_matsubara_measurements);
  } 
  else {
    itime_green_function_t bare_green_itime(n_tau+1, n_site, n_flavors);
    fourier_ptr_g0->backward_ft(bare_green_itime, bare_green_matsubara);
    evaluate_selfenergy_measurement_itime_rs(results, green_itime_measured, bare_green_itime, 
                                             beta, n_site, n_flavors, n_tau, n_self);
  }
  //Fourier transformations
  const auto pair_mean = results.at("n_i n_j").mean();
  const std::vector<double> density_pairs(pair_mean.data(), pair_mean.data()+pair_mean.size());
  FourierTransformer::generate_transformer_U(parms, fourier_ptr, densities,
                                            U_matrix(parms, input), density_pairs);
  if (measure_in_matsubara) {
    fourier_ptr->append_tail(green_matsubara_measured, bare_green_matsubara, n_matsubara_measurements);
    fourier_ptr->backward_ft(green_itime_measured, green_matsubara_measured);
  }
  else 
    fourier_ptr->forward_ft(green_itime_measured, green_matsubara_measured);
  require_finite(green_matsubara_measured, "CT-INT G_omega");
  require_finite(green_itime_measured, "CT-INT G_tau");
  green_matsubara_measured.write_hdf5(archive, "/G_omega");
  green_itime_measured.write_hdf5(archive, "/G_tau");
} 
