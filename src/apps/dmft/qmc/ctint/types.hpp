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

#include <alps/numeric/matrix.hpp>
#include "../types.h"
#include "operator.hpp"

/*types*/
class c_or_cdagger;
class vertex;
typedef class histogram simple_hist;
typedef std::vector<vertex> vertex_array;



enum measurement_methods {
  selfenergy_measurement_matsubara, //measurement using self energy method
  selfenergy_measurement_itime_rs, //measurement using self energy method in imag time, real space
};



class histogram
{
public:

  histogram(unsigned int N):hist_(N, 0){}

  unsigned long &operator[](unsigned int n){return hist_[n];}
  const unsigned long &operator[](unsigned int n) const{return hist_[n];}
  unsigned int size() const{return hist_.size();}

  unsigned int max_index() const
  {
    unsigned int max_index=0;
    double max=0;
    for(unsigned int i=0;i<hist_.size();++i){
      if(max<hist_[i]){
        max=hist_[i];
        max_index=i;
      }
    }
    return max_index;
  }

  unsigned int top_index() const
  {
    unsigned int top_index=0;
    for(unsigned int i=0;i<hist_.size();++i){
      if(hist_[i]!=0){
        top_index=i;
      }
    }
    return top_index;
  }

  double max(const unsigned int index) const
  {
    double max=0;
    for(unsigned int i=0;i<index;++i){
      if(max<hist_[i]){
        max=hist_[i];
      }
    }
    return max;
  }

  double average(const unsigned int index) const{
    double average=0;
    for(unsigned int i=0;i<index;++i){
      average+=hist_[i];
    }
    return average/index;
  }

  bool is_flat(const unsigned int index)const{return max(index)*0.8<average(index);}

  void clear()
  {
    for(unsigned int i=0;i<hist_.size();++i){
      hist_[i]=0;
    }
  }

private:

  std::vector<unsigned long> hist_;
};



typedef class vertex
{
public:
  vertex(const spin_t &flavor1, const site_t &site1, const unsigned int &c_dagger_1, const unsigned int &c_1,
         const spin_t &flavor2, const site_t &site2, const unsigned int &c_dagger_2, const unsigned int &c_2,
         const double &abs_w)
  {
    z1_=flavor1;
    z2_=flavor2;
    s1_=site1;
    s2_=site2;
    c1dagger_=c_dagger_1;
    c2dagger_=c_dagger_2;
    c1_=c_1;
    c2_=c_2;
    abs_w_=abs_w;
  }

  inline const double &abs_w() const {return abs_w_;}
  inline const unsigned int &flavor1() const {return z1_;}
  inline const unsigned int &flavor2() const {return z2_;}
  inline const unsigned int &site1() const {return s1_;}
  inline const unsigned int &site2() const {return s2_;}
  inline void set_site1(site_t site1) {s1_=site1;}
  inline void set_site2(site_t site2) {s2_=site2;}
  inline const unsigned int &c_dagger_1() const {return c1dagger_;}
  inline const unsigned int &c_dagger_2() const {return c2dagger_;}
  inline const unsigned int &c_1() const {return c1_;}
  inline const unsigned int &c_2() const {return c2_;}
  inline unsigned int &flavor1() {return z1_;}
  inline unsigned int &flavor2() {return z2_;}
  inline unsigned int &c_dagger_1() {return c1dagger_;}
  inline unsigned int &c_dagger_2() {return c2dagger_;}
  inline unsigned int &c_1() {return c1_;}
  inline unsigned int &c_2() {return c2_;}
private:
  unsigned int z1_, z2_;
  unsigned int s1_, s2_;
  unsigned int c1_, c2_;
  unsigned int c1dagger_, c2dagger_;
  double abs_w_;
} vertex;



class inverse_m_matrix
{
public:
  alps::numeric::matrix<double> &matrix() { return matrix_;}
  alps::numeric::matrix<double> const &matrix() const { return matrix_;}
  std::vector<creator> &creators(){ return creators_;}
  const std::vector<creator> &creators() const{ return creators_;}
  std::vector<annihilator> &annihilators(){ return annihilators_;}
  const std::vector<annihilator> &annihilators()const{ return annihilators_;}
  std::vector<double> &alpha(){ return alpha_;}
  const std::vector<double> &alpha() const{ return alpha_;}
private:
  alps::numeric::matrix<double> matrix_;
  std::vector<creator> creators_;         //an array of creation operators c_dagger corresponding to the row of the matrix
  std::vector<annihilator> annihilators_; //an array of to annihilation operators c corresponding to the column of the matrix
  std::vector<double> alpha_;             //an array of doubles corresponding to the alphas of Rubtsov for the c, cdaggers at the same index.
};
