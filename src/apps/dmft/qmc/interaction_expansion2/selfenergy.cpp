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
#include <complex>




void InteractionExpansion::compute_W_matsubara()
{
  std::vector<std::vector<std::valarray<std::complex<double> > > >Wk(n_flavors);
  for(unsigned int z=0;z<n_flavors;++z){
    Wk[z].resize(n_site);
    for(unsigned int j=0;j<n_site;++j){
      Wk[z][j].resize(n_matsubara);
      memset(&(Wk[z][j][0]), 0, sizeof(std::complex<double>)*(n_matsubara));
    }
  }
  measure_Wk(Wk, n_matsubara_measurements);
}



void InteractionExpansion::measure_Wk(std::vector<std::vector<std::valarray<std::complex<double> > > >& Wk,
                                         const unsigned int nfreq)
{
  for (unsigned int z=0; z<n_flavors; ++z) {
    assert( num_rows(M[z].matrix()) == num_cols(M[z].matrix()) );
    for (unsigned int k=0; k<n_site; k++) {
      for(unsigned int p=0;p<num_rows(M[z].matrix());++p){
        M[z].creators()[p].compute_exp(n_matsubara, +1);
        for(unsigned int q=0;q<num_cols(M[z].matrix());++q){
          M[z].annihilators()[q].compute_exp(n_matsubara, -1);
          std::complex<double>* Wk_z_k1_k2 = &Wk[z][k][0];
          const std::complex<double>* exparray_creators = M[z].creators()[p].exp_iomegat();
          const std::complex<double>* exparray_annihilators = M[z].annihilators()[q].exp_iomegat();
          std::complex<double> tmp = M[z].matrix()(p,q);
#pragma ivdep
          for(unsigned int o=0; o<nfreq; ++o){
            *Wk_z_k1_k2++ += (*exparray_creators++)*(*exparray_annihilators++)*tmp;
          }

        }
      }
    }
  }
  for(unsigned int flavor=0;flavor<n_flavors;++flavor){
    for (unsigned int k=0; k<n_site; k++) {
      std::stringstream Wk_real_name, Wk_imag_name;
      Wk_real_name  << "Wk_real_"  << flavor << "_" << k << "_" << k;
      Wk_imag_name  << "Wk_imag_"  << flavor << "_" << k << "_" << k;
      std::valarray<double> Wk_real(nfreq);
      std::valarray<double> Wk_imag(nfreq);
      for (unsigned int w=0; w<nfreq; ++w) {
        Wk_real[w] = Wk[flavor][k][w].real();
        Wk_imag[w] = Wk[flavor][k][w].imag();
      }
      record_measurement(Wk_real_name.str().c_str(), static_cast<std::valarray<double> > (Wk_real*sign));
      record_measurement(Wk_imag_name.str().c_str(), static_cast<std::valarray<double> > (Wk_imag*sign));
    }
  }
}



void InteractionExpansion::measure_densities()
{
  const double tau = beta*random();
  std::valarray<double> densities(n_flavors);
  for (spin_t flavor = 0; flavor < n_flavors; ++flavor) {
    const auto& matrix = M[flavor].matrix();
    alps::numeric::vector<double> g0_j(num_rows(matrix)), M_g0_j(num_rows(matrix));
    for (std::size_t j = 0; j < num_rows(matrix); ++j)
      g0_j[j] = green0_spline(M[flavor].creators()[j].t()-tau, flavor);
    if (num_rows(matrix)) gemv(matrix, g0_j, M_g0_j);
    double density = 1 + green0_spline(0, flavor);
    for (std::size_t i = 0; i < num_rows(matrix); ++i)
      density -= green0_spline(tau-M[flavor].annihilators()[i].t(), flavor) * M_g0_j[i];
    densities[flavor] = density;
    const auto name = (measurement_method == selfenergy_measurement_itime_rs ?
                      "density_" : "densities_") + std::to_string(flavor);
    record_measurement(name, density*sign);
  }
  record_measurement("densities", static_cast<std::valarray<double>>(densities*sign));
  std::valarray<double> ninj(std::size_t(n_flavors)*n_flavors);
  for (spin_t i = 0; i < n_flavors; ++i)
    for (spin_t j = 0; j < n_flavors; ++j)
      // Conditional estimates factorize between distinct flavors; n_i^2 = n_i.
      ninj[std::size_t(i)*n_flavors+j] = i == j ? densities[i] : densities[i]*densities[j];
  record_measurement("n_i n_j", static_cast<std::valarray<double>>(ninj*sign));
  record_measurement("density_correlation", densities[0]*densities[1]*sign);
  if (measurement_method == selfenergy_measurement_itime_rs && n_flavors == 2) {
    record_measurement("Sz_0", (densities[0]-densities[1])*sign);
    const double spin_squared = densities[0]+densities[1]-2*densities[0]*densities[1];
    record_measurement("Sz2_0", spin_squared*sign);
    record_measurement("Sz0_Sz0", spin_squared*sign);
  }
}



void InteractionExpansion::compute_W_itime()
{
  std::vector<std::vector<std::vector<std::valarray<double> > > >W_z_i_j(n_flavors);
  //first index: flavor. Second index: momentum. Third index: self energy tau point.
  for(unsigned int z=0;z<n_flavors;++z){
    W_z_i_j[z].resize(n_site);
    for(unsigned int i=0;i<n_site;++i){
      W_z_i_j[z][i].resize(n_site);
      for(unsigned int j=0;j<n_site;++j){
        W_z_i_j[z][i][j].resize(n_self+1);
        memset(&(W_z_i_j[z][i][j][0]), 0, sizeof(double)*(n_self+1));
      }
    }
  }
  int ntaupoints=10; //# of tau points at which we measure.
  std::vector<double> tau_2(ntaupoints);
  for(int i=0; i<ntaupoints;++i)
    tau_2[i]=beta*random();
  for(unsigned int z=0;z<n_flavors;++z){                  //loop over flavor
    assert( num_rows(M[z].matrix()) == num_cols(M[z].matrix()) );
    alps::numeric::matrix<double> g0_tauj(num_cols(M[z].matrix()),ntaupoints);
    alps::numeric::matrix<double> M_g0_tauj(num_rows(M[z].matrix()),ntaupoints);
    for(unsigned int s2=0;s2<n_site;++s2){             //site loop - second site.
      for(int j=0;j<ntaupoints;++j) {
        for(unsigned int i=0;i<num_cols(M[z].matrix());++i){ //G0_{s_p s_2}(tau_p - tau_2) where we set t2=0.
            g0_tauj(i,j) = green0_spline(M[z].creators()[i].t()-tau_2[j], z, M[z].creators()[i].s(), s2);
        }
      }
      if (num_rows(M[z].matrix())>0)
          gemm(M[z].matrix(), g0_tauj, M_g0_tauj);
      for(int j=0;j<ntaupoints;++j) {
        for(unsigned int p=0;p<num_rows(M[z].matrix());++p){       //operator one
          double sgn=1;
          double delta_tau=M[z].creators()[p].t()-tau_2[j];
          if(delta_tau<0){
            sgn=-1;
            delta_tau+=beta;
          }
          int bin=(int)(delta_tau/beta*n_self+0.5);
          site_t site_p=M[z].creators()[p].s();
          W_z_i_j[z][site_p][s2][bin] += M_g0_tauj(p,j)*sgn;
        }
      }
    }
  }
  if(is_thermalized()){
    for(unsigned int flavor=0;flavor<n_flavors;++flavor){
      for(unsigned int i=0;i<n_site;++i){
        for(unsigned int j=0;j<n_site;++j){
          std::stringstream W_name;
          W_name  <<"W_"  <<flavor<<"_"<<i<<"_"<<j;
          record_measurement(W_name  .str().c_str(), static_cast<std::valarray<double> > (W_z_i_j[flavor][i][j]*(sign/ntaupoints)));
        }
      }
    }
  }
}



void evaluate_selfenergy_measurement_matsubara(const InteractionExpansion::results_type &results,
                                                                        matsubara_green_function_t &green_matsubara_measured,
                                                                        const matsubara_green_function_t &bare_green_matsubara,
                                                                        const double &beta, std::size_t n_site,
                                                                        std::size_t n_flavors, std::size_t n_matsubara)
{
  double max_error = 0.;
  std::cout<<"evaluating self energy measurement: matsubara, reciprocal space"<<std::endl;
  matsubara_green_function_t Wk(n_matsubara, n_site, n_flavors);
  Wk.clear();
  for(std::size_t z=0;z<n_flavors;++z){
    for (std::size_t k=0; k<n_site; k++) {
      std::stringstream Wk_real_name, Wk_imag_name;
      Wk_real_name  <<"Wk_real_"  <<z<<"_"<<k << "_" << k;
      Wk_imag_name  <<"Wk_imag_"  <<z<<"_"<<k << "_" << k;
      auto mean_real = results.at(Wk_real_name.str().c_str()).mean();
      auto mean_imag = results.at(Wk_imag_name.str().c_str()).mean();
      for(unsigned int w=0;w<n_matsubara;++w)
        Wk(w, k, k, z) = std::complex<double>(mean_real[w], mean_imag[w])/( beta*n_site);
      auto error_real = results.at(Wk_real_name.str().c_str()).stderror();
      auto error_imag = results.at(Wk_imag_name.str().c_str()).stderror();
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
  for(std::size_t z=0;z<n_flavors;++z)
    for (std::size_t k=0; k<n_site; k++)
      for(std::size_t w=0;w<n_matsubara;++w)
        green_matsubara_measured(w,k,k, z) = bare_green_matsubara(w,k,k,z)
        - bare_green_matsubara(w,k,k,z) * bare_green_matsubara(w,k,k,z) * Wk(w,k,k,z);
}


double green0_spline(const itime_green_function_t &green0, const itime_t delta_t,
                                              const int s1, const int s2, const spin_t z, int n_tau, double beta);


void evaluate_selfenergy_measurement_itime_rs(const InteractionExpansion::results_type &results,
                                                                       itime_green_function_t &green_result,
                                                                       const itime_green_function_t &green0,
                                                                       const double &beta, const int n_site,
                                                                       const int n_flavors, const int n_tau, const int n_self)
{
  std::cout<<"evaluating self energy measurement: itime, real space."<<std::endl;
  clock_t time0=clock();
  std::vector<std::vector<std::vector<std::vector<double> > > >W_z_i_j(n_flavors);
  //first index: flavor. Second index: momentum. Third index: self energy tau point.
  double max_error = 0.;
  for(int z=0;z<n_flavors;++z){
    W_z_i_j[z].resize(n_site);
    for(int i=0;i<n_site;++i){
      W_z_i_j[z][i].resize(n_site);
    }
    for(int i=0;i<n_site;++i){
      for(int j=0;j<n_site;++j){
        std::stringstream W_name;
        W_name<<"W_"<<z<<"_"<<i<<"_"<<j;
        W_z_i_j[z][i][j].resize(n_self+1);
        auto tmp = results.at(W_name.  str().c_str()).mean();
        auto errorvec = results.at(W_name.  str().c_str()).stderror();
        for(int k=0;k<n_self+1;++k){
          W_z_i_j[z][i][j][k]=tmp[k];
          double error = errorvec[k];
          max_error = error>max_error ? error : max_error;
        }
      }
    }
    std::cout << "Maximum error in W: " << max_error << std::endl;
    for(int i=0;i<n_site;++i){
      for(int j=0;j<n_site;++j){
        for(int k=0;k<n_tau;++k){
          green_result(k,i,j,z)=green0(k,i,j,z);
          double timek=(k/(double)n_tau)*beta;
          for(int l=0;l<=n_self;++l){
            double timel;
            if(l==0){
              timel=(0.5/(double)n_self)*beta; //middle position of our first bin.
            }else if(l==n_self){
              timel=((n_self-0.5)/(double)n_self)*beta; //middle position of our last bin.
            }else{
              timel=(l/(double)n_self)*beta; //middle position of our remaining bins.
            }
            for(int p=0;p<n_site;++p){
              //this will make a tiny error if the bins are equal.
              green_result(k,i,j,z)-=green0_spline(green0, timek-timel,i,p,z, n_tau, beta)*W_z_i_j[z][p][j][l];
            }
          }
        }
        green_result(n_tau,i,j,z)=(i==j?-1:0)-green_result(0,i,j,z);
      }
    }
  }
  clock_t time1=clock();
  std::cout<<"evaluate of SE measurement took: "<<(time1-time0)/(double)CLOCKS_PER_SEC<<std::endl;
}



template<class X, class Y> inline Y linear_interpolate(const X x0, const X x1, const Y y0, const Y y1, const X x)
{
  return y0 + (x-x0)/(x1-x0)*(y1-y0);
}





double green0_spline(const itime_green_function_t &green0, const itime_t delta_t,
                                              const int s1, const int s2, const spin_t z, int n_tau, double beta)
{
  double temperature=1./beta;
  double almost_zero=1.e-12;
  double n_tau_inv=1./n_tau;
  if(delta_t*delta_t < almost_zero){
    return green0(0,s1,s2,z);
  } else if(delta_t>0){
    int time_index_1 = (int)(delta_t*n_tau*temperature);
    int time_index_2 = time_index_1+1;
    return linear_interpolate((double)time_index_1*beta*n_tau_inv, (double)time_index_2*beta*n_tau_inv,green0(time_index_1,s1,s2,z),
                              green0(time_index_2,s1,s2,z),delta_t);
  } else{
    int time_index_1 = (int)((beta+delta_t)*n_tau*temperature);
    int time_index_2 = time_index_1+1;
    return -linear_interpolate((double)time_index_1*beta*n_tau_inv, (double)time_index_2*beta*n_tau_inv, green0(time_index_1,s1,s2,z),
                               green0(time_index_2, s1,s2,z), delta_t+beta);
  }
}
