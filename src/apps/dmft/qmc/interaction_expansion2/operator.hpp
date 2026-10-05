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

#ifndef DMFT_QMC_WEAK_COUPLING_OPERATOR_H
#define DMFT_QMC_WEAK_COUPLING_OPERATOR_H



/*creation and annihilation operator class*/
typedef class c_or_cdagger   //represents a creation operator or an annihilation operator
{ 
public:
  c_or_cdagger( const spin_t z,const site_t s, const itime_t t, const frequency_t n_matsubara)
  {
    s_ = s;
    z_ = z;
    t_ = t;
    nm_ = n_matsubara;
    exp_computed_ = false;
    exp_iomegat_ = 0;
  }
  
  
  
  ~c_or_cdagger()
  {
    if(exp_computed_)
      delete [] exp_iomegat_;
  }
  
  
  
  const c_or_cdagger & operator=(const c_or_cdagger &c)
  {
    if(this != &c){
      s_ = c.s_;
      z_ = c.z_;
      t_ = c.t_;
      if(use_static_exp_)
        exp_iomegat_=c.exp_iomegat_;
      else {
         if (exp_computed_ && c.exp_computed_) {
          memcpy(exp_iomegat_, c.exp_iomegat_, sizeof(std::complex<double>)*c.nm_);
        }
        else if (exp_computed_ && (!c.exp_computed_)) {
          delete [] exp_iomegat_;
        }
        else if ((!exp_computed_) && c.exp_computed_) {
          exp_iomegat_ = new std::complex<double>[c.nm_];
          memcpy(exp_iomegat_, c.exp_iomegat_, sizeof(std::complex<double>)*c.nm_);
        }
      }
      nm_=c.nm_;
      exp_computed_=c.exp_computed_;          
    }
    return *this;
  }
  
  
  
  c_or_cdagger(const c_or_cdagger & c)
  {
    exp_computed_=false;
    operator=(c);
  }
  


  inline const spin_t &flavor() const {return z_;}
  inline spin_t &flavor() {return z_;}
  inline const itime_t &t() const{return t_;}
  inline itime_t &t() {return t_;}
  inline const site_t &s() const {return s_;}
  inline void flavor(spin_t z){z_=z;}
  inline void s(site_t s){s_=s;}
  inline const std::complex<double> * exp_iomegat() const {return exp_iomegat_;} 
  //contains exp(iomegat) if its a creator, exp(-iomegat) if it's an annihilator.
  static void initialize_simulation(const alps::params &parms);
  
  
  static const std::complex<double> *exp_iomegan_tau(const double &tau) 
  {
    int taun=(int)(tau*ntau_/beta_); 
    return &(exp_iomegan_tau_[std::size_t(taun)*2*nm_]);
  }
  

  static const std::complex<double> *exp_min_iomegan_tau(const double &tau) 
  {
    int taun=(int)(tau*ntau_/beta_);
    return &(exp_iomegan_tau_[std::size_t(taun)*2*nm_ + nm_]);
  }
  
  
  
private:
  site_t s_;      //this vertex's site
  itime_t t_;     //its imaginary time point
  spin_t z_;      //its flavor or flavor
  static unsigned int nm_;        //number of matsubara frequencies
  std::complex<double> *exp_iomegat_;
  bool exp_computed_;
  static bool use_static_exp_; //do we compute the exps once only or for each time slice?
  static unsigned int ntau_;
  static double beta_;
  static std::vector<double> omegan_;
  static std::vector<std::complex<double>> exp_iomegan_tau_;



public:
  

  void compute_exp(const frequency_t n_matsubara, const int sign)
  {
    if(!use_static_exp_){
      if(!exp_computed_){
        exp_iomegat_=new std::complex<double>[n_matsubara];
        for(frequency_t o=0;o<n_matsubara;++o) {
          const double phase=omegan_[o]*t_;
          exp_iomegat_[o]={std::cos(phase), sign*std::sin(phase)};
        }
        exp_computed_=true;
      }
    } else { //use static exp
      int taun=(int)(t_*ntau_/beta_);
      assert(taun<ntau_);
      if(sign==1) 
        exp_iomegat_=&(exp_iomegan_tau_[std::size_t(taun)*2*nm_]);
      else
        exp_iomegat_=&(exp_iomegan_tau_[std::size_t(taun)*2*nm_ + nm_]);
    }


  }

} creator, annihilator;

#endif
