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
#include <alps/alea.h>
#include "alps/ngs/make_deprecated_parameters.hpp"

void evaluate_selfenergy_measurement_matsubara(const alps::results_type<HubbardInteractionExpansion>::type &results, 
                                                                        matsubara_green_function_t &green_matsubara_measured,
                                                                        const matsubara_green_function_t &bare_green_matsubara, 
                                                                        std::vector<double>& densities,
                                                                        const double &beta, std::size_t n_site, 
                                                                        std::size_t n_flavors, std::size_t n_matsubara);
void evaluate_selfenergy_measurement_itime_rs(const alps::results_type<HubbardInteractionExpansion>::type &results, 
                                                                       itime_green_function_t &green_result,
                                                                       const itime_green_function_t &green0, 
                                                                       const double &beta, const int n_site, 
                                                                       const int n_flavors, const int n_tau, const int n_self);



void compute_greens_functions(const alps::results_type<HubbardInteractionExpansion>::type &results, const alps::parameters_type<HubbardInteractionExpansion>::type& parms, const alps::params &input, const alps::params &output)
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
  FourierTransformer::generate_transformer(alps::make_deprecated_parameters(parms), fourier_ptr_g0);
  //find whether our data is in imaginary time or frequency:
  bool measure_in_matsubara=true;
  if(parms.value_or("HISTOGRAM_MEASUREMENT", false))
    measure_in_matsubara=false;
  std::vector<double> mean_order=results["PertOrder"].mean<std::vector<double> >();
  
  std::cout<<"average matrix size was: "<<std::endl;
  std::ofstream matrix_size;
  if (output.exists("matrix_size")) {
    matrix_size.open(output["matrix_size"].as<std::string>(), std::ios::app);
    if (!matrix_size)
      throw std::runtime_error("cannot open CT-INT matrix_size output");
  }
  for(unsigned int i=0;i<n_flavors;++i){
    std::cout<<mean_order[i]<<"\t";
    if (matrix_size.is_open()) matrix_size<<mean_order[i]<<"\t";
  }
  std::cout<<std::endl;
  if (matrix_size.is_open()) matrix_size<<std::endl;
  std::cout<<"average sign was: "<<results["Sign"].mean<double>()<<" error: "<<results["Sign"].error<double>()<<std::endl;
  //single particle Green function measurements
  matsubara_green_function_t bare_green_matsubara(n_matsubara, n_site, n_flavors);
  std::vector<double> densities(n_flavors);
  read_ctint_bare_green(parms, input, bare_green_matsubara);
  if(measure_in_matsubara) {
    evaluate_selfenergy_measurement_matsubara(results, green_matsubara_measured, 
                                              bare_green_matsubara, densities, 
                                              beta, n_site, n_flavors, n_matsubara_measurements);
  } 
  else {
    itime_green_function_t bare_green_itime(n_tau+1, n_site, n_flavors);
    fourier_ptr_g0->backward_ft(bare_green_itime, bare_green_matsubara);
    evaluate_selfenergy_measurement_itime_rs(results, green_itime_measured, bare_green_itime, 
                                             beta, n_site, n_flavors, n_tau, n_self);
  }
  //Fourier transformations
  if (!measure_in_matsubara) {
    for (unsigned int z=0; z<n_flavors; ++z) {
      densities[z] = 0;
      for (unsigned int i=0; i<n_site; ++i)
        densities[z] -= green_itime_measured(n_tau,i,i,z);
      densities[z] /= n_site;
    }
  }
  FourierTransformer::generate_transformer_U(alps::make_deprecated_parameters(parms), fourier_ptr, densities);
  if (measure_in_matsubara) {
    fourier_ptr->append_tail(green_matsubara_measured, bare_green_matsubara, n_matsubara_measurements);
    fourier_ptr->backward_ft(green_itime_measured, green_matsubara_measured);
  }
  else 
    fourier_ptr->forward_ft(green_itime_measured, green_matsubara_measured);
  {
    alps::hdf5::archive ar(output["results"].as<std::string>(), "a");
    green_matsubara_measured.write_hdf5(ar, "/G_omega");
    green_itime_measured.write_hdf5(ar, "/G_tau");
  }
} 
