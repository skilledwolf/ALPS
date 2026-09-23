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

#pragma once
#include "update_kernel.hpp"
#include "../green_function.h"
#include <complex>
#include <cstring>
#include <sstream>
#include <valarray>
#ifdef SSE
#include <emmintrin.h>
#endif

namespace alps { namespace ctint {

#ifdef SSE
class twocomplex{
public:
  inline twocomplex(){};
  inline twocomplex(const std::complex<double>&p, const std::complex<double> &q){
    r=_mm_loadl_pd(r, &(p.real())); //load two complex numbers.
    r=_mm_loadh_pd(r, &(q.real()));
    i=_mm_loadl_pd(i, &(p.imag()));
    i=_mm_loadh_pd(i, &(q.imag()));
  }
  inline void store(std::complex<double> &p, std::complex<double> &q){
    _mm_store_sd(&(p.real()), r);
    _mm_store_sd(&(p.imag()), i);
    _mm_storeh_pd(&(q.real()), r);
    _mm_storeh_pd(&(q.imag()), i);
  }
  __m128d r;
  __m128d i;
};

inline twocomplex fastcmult(const twocomplex &a, const twocomplex &b)
{
  twocomplex c;
  c.r = _mm_sub_pd(_mm_mul_pd(a.r,b.r), _mm_mul_pd(a.i,b.i));
  c.i = _mm_add_pd(_mm_mul_pd(a.r,b.i), _mm_mul_pd(a.i,b.r));
  return c;
}
#endif

// Runtime adapters own observable types and names; these kernels retain the
// operation and random-draw order of both solvers.
struct measurement_kernel {
  template<class Run>
  static void compute_W_matsubara(Run& run)
  {
    static std::vector<std::vector<std::valarray<std::complex<double> > > >Wk(run.n_flavors);
    for(unsigned int z=0;z<run.n_flavors;++z){
      Wk[z].resize(run.n_site);
      for(unsigned int j=0;j<run.n_site;++j){
        Wk[z][j].resize(run.n_matsubara);
        memset(&(Wk[z][j][0]), 0, sizeof(std::complex<double>)*(run.n_matsubara));
      }
    }
    measure_Wk(run, Wk, run.n_matsubara_measurements);
    measure_densities(run);
  }

  template<class Run>
  static void measure_Wk(Run& run, std::vector<std::vector<std::valarray<std::complex<double> > > >& Wk,
                                           const unsigned int nfreq)
  {
    for (unsigned int z=0; z<run.n_flavors; ++z) {
      assert( num_rows(run.M[z].matrix()) == num_cols(run.M[z].matrix()) );
      for (unsigned int k=0; k<run.n_site; k++) {
        for(unsigned int p=0;p<num_rows(run.M[z].matrix());++p){
          run.M[z].creators()[p].compute_exp(run.n_matsubara, +1);
          for(unsigned int q=0;q<num_cols(run.M[z].matrix());++q){
            run.M[z].annihilators()[q].compute_exp(run.n_matsubara, -1);
            std::complex<double>* Wk_z_k1_k2 = &Wk[z][k][0];
            const std::complex<double>* exparray_creators = run.M[z].creators()[p].exp_iomegat();
            const std::complex<double>* exparray_annihilators = run.M[z].annihilators()[q].exp_iomegat();
            std::complex<double> tmp = run.M[z].matrix()(p,q);
#ifndef SSE
  //#pragma ivdep
            for(unsigned int o=0; o<nfreq; ++o){
              *Wk_z_k1_k2++ += (*exparray_creators++)*(*exparray_annihilators++)*tmp;
            }
#else
  //#pragma ivdep
            for(int o=0;o<nfreq;o+=2) {
              twocomplex exp_c(*exparray_creators++,    *exparray_creators++); //load it all into xmm registers
              twocomplex exp_a(*exparray_annihilators++,*exparray_annihilators++);
              twocomplex tmp2(tmp,tmp);
              twocomplex product=fastcmult(fastcmult(exp_c,exp_a),tmp2);
              std::complex<double> Wk1, Wk2;
              product.store(Wk1, Wk2);
              *Wk_z_k1_k2++ += Wk1;
              *Wk_z_k1_k2++ += Wk2;
            }
#endif
          }
        }
      }
    }
    for(unsigned int flavor=0;flavor<run.n_flavors;++flavor){
      for (unsigned int k=0; k<run.n_site; k++) {
        std::stringstream Wk_real_name, Wk_imag_name;
        Wk_real_name  << "Wk_real_"  << flavor << "_" << k << "_" << k;
        Wk_imag_name  << "Wk_imag_"  << flavor << "_" << k << "_" << k;
        std::valarray<double> Wk_real(nfreq);
        std::valarray<double> Wk_imag(nfreq);
        for (unsigned int w=0; w<nfreq; ++w) {
          Wk_real[w] = Wk[flavor][k][w].real();
          Wk_imag[w] = Wk[flavor][k][w].imag();
        }
        run.ctint_measure(Wk_real_name.str().c_str(), Wk_real*run.sign);
        run.ctint_measure(Wk_imag_name.str().c_str(), Wk_imag*run.sign);
      }
    }
  }

  template<class Run>
  static void measure_densities(Run& run)
  {
    std::vector< std::vector<double> > dens(run.n_flavors);
    for(unsigned int z=0;z<run.n_flavors;++z){
      dens[z].resize(run.n_site);
      memset(&(dens[z][0]), 0., sizeof(double)*(run.n_site));
    }
    double tau = run.beta*run.ctint_uniform();
    for (unsigned int z=0; z<run.n_flavors; ++z) {
      alps::numeric::vector<double> g0_tauj(num_rows(run.M[z].matrix()));
      alps::numeric::vector<double> M_g0_tauj(num_rows(run.M[z].matrix()));
      alps::numeric::vector<double> g0_taui(num_rows(run.M[z].matrix()));
      for (unsigned int s=0;s<run.n_site;++s) {
        for (unsigned int j=0;j<num_rows(run.M[z].matrix());++j)
          g0_tauj[j] = run.green0_spline(run.M[z].creators()[j].t()-tau, z, run.M[z].creators()[j].s(), s);
        for (unsigned int i=0;i<num_rows(run.M[z].matrix());++i)
          g0_taui[i] = run.green0_spline(tau-run.M[z].annihilators()[i].t(),z, s, run.M[z].annihilators()[i].s());
        if (num_rows(run.M[z].matrix())>0)
            gemv(run.M[z].matrix(),g0_tauj,M_g0_tauj);
        dens[z][s] += run.green0_spline(0,z,s,s);
        for (unsigned int j=0;j<num_rows(run.M[z].matrix());++j)
          dens[z][s] -= g0_taui[j]*M_g0_tauj[j];
      }
    }
    std::valarray<double> densities(0., run.n_flavors);
    for (unsigned int z=0; z<run.n_flavors; ++z) {
      std::valarray<double> densmeas(run.n_site);
      for (unsigned int i=0; i<run.n_site; ++i) {
        densities[z] += dens[z][i];
        densmeas[i] = 1+dens[z][i];
      }
      run.ctint_measure("densities_"+boost::lexical_cast<std::string>(z), densmeas*run.sign);
      densities[z] /= run.n_site;
      densities[z] = 1 + densities[z];
    }
    run.ctint_measure("densities", densities*run.sign);
    double density_correlation = 0.;
    for (unsigned int i=0; i<run.n_site; ++i) {
      density_correlation += (1+dens[0][i])*(1+dens[1][i]);
    }
    density_correlation /= run.n_site;
    run.measurements["density_correlation"] << density_correlation*run.sign;
    std::valarray<double> ninj(run.n_site*run.n_site*4);
    for (unsigned int i=0; i<run.n_site; ++i) {
      for (unsigned int j=0; j<run.n_site; ++j) {
        ninj[i*run.n_site+j] = (1+dens[0][i])*(1+dens[0][j]);
        ninj[i*run.n_site+j+1] = (1+dens[0][i])*(1+dens[1][j]);
        ninj[i*run.n_site+j+2] = (1+dens[1][i])*(1+dens[0][j]);
        ninj[i*run.n_site+j+3] = (1+dens[1][i])*(1+dens[1][j]);
      }
    }
    run.ctint_measure("n_i n_j", ninj*run.sign);
  }

  template<class Run>
  static void compute_W_itime(Run& run)
  {
    static std::vector<std::vector<std::vector<std::valarray<double> > > >W_z_i_j(run.n_flavors);
    //first index: flavor. Second index: momentum. Third index: self energy tau point.
    std::vector<std::vector<double> >density(run.n_flavors);
    for(unsigned int z=0;z<run.n_flavors;++z){
      W_z_i_j[z].resize(run.n_site);
      density[z].resize(run.n_site);
      for(unsigned int i=0;i<run.n_site;++i){
        W_z_i_j[z][i].resize(run.n_site);
        for(unsigned int j=0;j<run.n_site;++j){
          W_z_i_j[z][i][j].resize(run.n_self+1);
          memset(&(W_z_i_j[z][i][j][0]), 0, sizeof(double)*(run.n_self+1));
        }
      }
    }
    int ntaupoints=10; //# of tau points at which we measure.
    std::vector<double> tau_2(ntaupoints);
    for(int i=0; i<ntaupoints;++i)
      tau_2[i]=run.beta*run.ctint_uniform();
    for(unsigned int z=0;z<run.n_flavors;++z){                  //loop over flavor
      assert( num_rows(run.M[z].matrix()) == num_cols(run.M[z].matrix()) );
      alps::numeric::matrix<double> g0_tauj(num_cols(run.M[z].matrix()),ntaupoints);
      alps::numeric::matrix<double> M_g0_tauj(num_rows(run.M[z].matrix()),ntaupoints);
      alps::numeric::vector<double> g0_taui(num_rows(run.M[z].matrix()));
      for(unsigned int s2=0;s2<run.n_site;++s2){             //site loop - second site.
        for(int j=0;j<ntaupoints;++j) {
          for(unsigned int i=0;i<num_cols(run.M[z].matrix());++i){ //G0_{s_p s_2}(tau_p - tau_2) where we set t2=0.
              g0_tauj(i,j) = run.green0_spline(run.M[z].creators()[i].t()-tau_2[j], z, run.M[z].creators()[i].s(), s2);
          }
        }
        if (num_rows(run.M[z].matrix())>0)
            gemm(run.M[z].matrix(), g0_tauj, M_g0_tauj);
        for(int j=0;j<ntaupoints;++j) {
          for(unsigned int p=0;p<num_rows(run.M[z].matrix());++p){       //operator one
            double sgn=1;
            double delta_tau=run.M[z].creators()[p].t()-tau_2[j];
            if(delta_tau<0){
              sgn=-1;
              delta_tau+=run.beta;
            }
            int bin=(int)(delta_tau/run.beta*run.n_self+0.5);
            site_t site_p=run.M[z].creators()[p].s();
            W_z_i_j[z][site_p][s2][bin] += M_g0_tauj(p,j)*sgn;
          }
        }
        for(unsigned int i=0;i<num_rows(run.M[z].matrix());++i){
          g0_taui[i]=run.green0_spline(tau_2[0]-run.M[z].annihilators()[i].t(),z, s2, run.M[z].annihilators()[i].s());
        }
        density[z][s2]=run.green0_spline(0,z,s2,s2);
        for(unsigned int i=0;i<num_rows(run.M[z].matrix());++i){
          density[z][s2]-= g0_taui[i]*M_g0_tauj(i,0);
        }
      }
    }
    if(run.is_thermalized()){
      for(unsigned int flavor=0;flavor<run.n_flavors;++flavor){
        for(unsigned int i=0;i<run.n_site;++i){
          for(unsigned int j=0;j<run.n_site;++j){
            std::stringstream W_name;
            W_name  <<"W_"  <<flavor<<"_"<<i<<"_"<<j;
            run.ctint_measure(W_name  .str().c_str(), W_z_i_j[flavor][i][j]*(run.sign/ntaupoints));
          }
          std::stringstream density_name;
          density_name<<"density_"<<flavor;
          if (!Run::ctint_single_site_density || run.n_site > 1) density_name<<"_"<<i;
          run.measurements[density_name.str().c_str()]<<density[flavor][i]*run.sign;
          if(run.n_flavors==2 && flavor==0){ //then we know how to compute Sz^2
            std::stringstream sz_name, sz2_name, sz0_szj_name;
            sz_name<<"Sz_"<<i; sz2_name<<"Sz2_"<<i; sz0_szj_name<<"Sz0_Sz"<<i;
            run.measurements[sz_name.str().c_str()]<<(density[0][i]-density[1][i])*run.sign;
            run.measurements[sz2_name.str().c_str()]<<(density[0][i]-density[1][i])*(density[0][i]-density[1][i])*run.sign;
            run.measurements[sz0_szj_name.str().c_str()]<<(density[0][0]-density[1][0])*(density[0][i]-density[1][i])*run.sign;
          }
        }
      }
    }
  }

};

inline double green0_spline(const itime_green_function_t &green0, const itime_t delta_t,
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

template<class ReadW>
void evaluate_itime_green(itime_green_function_t& green_result,
    const itime_green_function_t& green0, double beta, int n_site,
    int n_flavors, int n_tau, int n_self, ReadW read_W) {
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
        auto measurement = read_W(W_name.str());
        auto const& tmp = measurement.first;
        auto const& errorvec = measurement.second;
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
}

} }
