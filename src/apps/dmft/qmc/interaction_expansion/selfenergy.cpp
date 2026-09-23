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
#include "../ctint/selfenergy.hpp"
#include <alps/alea.h>

void InteractionExpansionRun::ctint_measure(std::string const& name, std::valarray<double> const& value) {
  measurements.get<alps::SignedObservable<alps::RealVectorObservable> >(name) << value;
}

void InteractionExpansionRun::compute_W_matsubara() {
  alps::ctint::measurement_kernel::compute_W_matsubara(*this);
}

void InteractionExpansionRun::measure_Wk(std::vector<std::vector<std::valarray<std::complex<double> > > >& Wk, unsigned int nfreq) {
  alps::ctint::measurement_kernel::measure_Wk(*this, Wk, nfreq);
}

void InteractionExpansionRun::measure_densities() {
  alps::ctint::measurement_kernel::measure_densities(*this);
}

void InteractionExpansionRun::compute_W_itime() {
  alps::ctint::measurement_kernel::compute_W_itime(*this);
}

void InteractionExpansionSim::evaluate_selfenergy_measurement_matsubara(const alps::ObservableSet &gathered_measurements, 
                                                                        matsubara_green_function_t &green_matsubara_measured,
                                                                        const matsubara_green_function_t &bare_green_matsubara, 
                                                                        std::vector<double>& densities,
                                                                        const double &beta, const int n_site, 
                                                                        const int n_flavors, const int n_matsubara) const
{
  double max_error = 0.;
  std::cout<<"evaluating self energy measurement."<<std::endl;
  matsubara_green_function_t Wk(n_matsubara, n_site, n_flavors);
  Wk.clear();
  matsubara_green_function_t reduced_bare_green_matsubara(n_matsubara, n_site, n_flavors);
  reduced_bare_green_matsubara.clear();
  for(int z=0;z<n_flavors;++z){
    for (int k=0; k<n_site; k++) {                   
      std::stringstream Wk_real_name, Wk_imag_name;
      Wk_real_name  <<"Wk_real_"  <<z<<"_"<<k << "_" << k;
      Wk_imag_name  <<"Wk_imag_"  <<z<<"_"<<k << "_" << k;
      alps::RealVectorObsevaluator Weval_real=gathered_measurements[Wk_real_name.str().c_str()];
      alps::RealVectorObsevaluator Weval_imag=gathered_measurements[Wk_imag_name.str().c_str()];
      Weval_real /= beta*n_site;
      Weval_imag /= beta*n_site;
      std::valarray<double> mean_real = Weval_real.mean();
      std::valarray<double> mean_imag = Weval_imag.mean();
      for(int w=0;w<n_matsubara;++w)
        Wk(w, k, k, z) = std::complex<double>(mean_real[w], mean_imag[w]);
      for(int w=0;w<n_matsubara;++w)
        reduced_bare_green_matsubara(w, k, k, z) = bare_green_matsubara(w, k, k, z);
      std::valarray<double> error_real = Weval_real.error();
      std::valarray<double> error_imag = Weval_imag.error();
      for (unsigned int e=0; e<error_real.size(); ++e) {
        double ereal = error_real[e];
        double eimag = error_imag[e];
        double error = (ereal >= eimag) ? ereal : eimag;
        max_error = (error > max_error) ? error : max_error;
      }
    }
  }
  std::cout << "Maximal error in Wk: " << max_error << std::endl;
  green_matsubara_measured.clear();
  for(int z=0;z<n_flavors;++z)
    for (int k=0; k<n_site; k++)                    
      for(int w=0;w<n_matsubara;++w)
        green_matsubara_measured(w,k,k, z) = bare_green_matsubara(w,k,k,z) 
        - bare_green_matsubara(w,k,k,z) * bare_green_matsubara(w,k,k,z) * Wk(w,k,k,z);
  std::valarray<double> dens = alps::RealVectorObsevaluator(gathered_measurements["densities"]).mean();
  for (int z=0; z<n_flavors; ++z) 
    densities[z] = dens[z];
}




void InteractionExpansionSim::evaluate_selfenergy_measurement_itime_rs(const alps::ObservableSet &gathered_measurements, 
                                                                       itime_green_function_t &green_result,
                                                                       const itime_green_function_t &green0, 
                                                                       const double &beta, const int n_site, 
                                                                       const int n_flavors, const int n_tau, const int n_self) const
{
  std::cout<<"evaluating self energy measurement: itime, real space."<<std::endl;
  clock_t time0=clock();
  alps::ctint::evaluate_itime_green(green_result, green0, beta, n_site,
    n_flavors, n_tau, n_self, [&](std::string const& name) {
      alps::RealVectorObsevaluator value = gathered_measurements[name];
      return std::make_pair(value.mean(), value.error());
    });
  if(n_flavors==2){  
    alps::RealObsevaluator Sz_obs = gathered_measurements["Sz_0"];
    for(int i=1;i<n_site;i++) {
      std::stringstream Sz_name;
      Sz_name << "Sz_" << i;
      alps::RealObsevaluator Sz_i_obs = gathered_measurements[Sz_name.str().c_str()];
      Sz_i_obs *= ((i%2==0) ? 1 : -1);
      Sz_obs += Sz_i_obs;
    }
    Sz_obs /= n_site;
    std::ofstream szstream("staggered_sz", std::ios::app);
    szstream << Sz_obs.mean() << "\t" << Sz_obs.error() << std::endl;
  }
  clock_t time1=clock();
  std::cout<<"evaluate of SE measurement took: "<<(time1-time0)/(double)CLOCKS_PER_SEC<<std::endl;
}



double InteractionExpansionSim::green0_spline(const itime_green_function_t &green0, const itime_t delta_t, 
                                              const int s1, const int s2, const spin_t z, int n_tau, double beta) const
{
  return alps::ctint::green0_spline(green0, delta_t, s1, s2, z, n_tau, beta);
}
