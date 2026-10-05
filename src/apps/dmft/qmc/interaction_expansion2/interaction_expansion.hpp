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

#ifndef DMFT_QMC_WEAK_COUPLING_H
#define DMFT_QMC_WEAK_COUPLING_H

#include <alps/alea/batch.hpp>
#include <alps/ngs/random01.hpp>
#include <alps/run_config.hpp>

#include <functional>
#include <map>
#include <cmath>
#include "green_function.h"
#include "types.h"
#include "fouriertransform.h"
#include "U_matrix.h"
#include "operator.hpp"
#include <alps/numeric/matrix.hpp>


/*types*/
class c_or_cdagger;
class vertex;
typedef std::vector<vertex> vertex_array;



enum measurement_methods {
  selfenergy_measurement_matsubara, //measurement using self energy method
  selfenergy_measurement_itime_rs, //measurement using self energy method in imag time, real space
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

class InteractionExpansion
{
public:

  InteractionExpansion(const alps::run_configuration& run, int rank);
  using results_type = std::map<std::string, alps::alea::batch_result<double>>;
  bool run(std::function<bool()> const& stop_callback);
  results_type collect_results(alps::alea::reducer const* reduction = nullptr) const;
  bool is_thermalized() const {return step > therm_steps;}
  void update();
  void measure();
  double fraction_completed() const;

protected:

  struct measurement {
    bool signed_value;
    alps::alea::batch_acc<double> accumulator;
  };
  alps::random01 random;
  std::map<std::string, measurement> measurements;
  void record_measurement(std::string const&, std::valarray<double> const&);
  void record_measurement(std::string const&, double);

  /*functions*/
  /*green's function*/
  // in file spines.cpp
  double green0_spline(const c_or_cdagger &cdagger, const c_or_cdagger &c) const;
  double green0_spline(const itime_t delta_t, const spin_t flavor, const site_t site1, const site_t site2) const;
  double green0_spline(const itime_t delta_t, const spin_t flavor) const;

  /*the actual solver functions*/
  // in file solver.cpp
  void interaction_expansion_step(void);
  void reset_perturbation_series(void);

  // in file fastupdate.cpp:
  double fastupdate_up(const int operator_nr, bool compute_only_weight);
  double fastupdate_down(const int operator_nr, const int flavor, bool compute_only_weight);

  /*measurement functions*/
  // in file measurements.cpp
  void measure_observables(void);
  void initialize_observables(void);

  void compute_W_matsubara();
  void compute_W_itime();
  void measure_Wk(std::vector<std::vector<std::valarray<std::complex<double> > > >& Wk, const unsigned int nfreq);
  void measure_densities();

  double try_add();
  void perform_add();
  void reject_add();
  double try_remove(unsigned int vertex_nr);
  void perform_remove(unsigned int vertex_nr);

  const std::size_t num_bins;

  /*private member variables, constant throughout the simulation*/
  const unsigned int max_order;
  const spin_t n_flavors;                                //number of flavors (called 'flavors') in InteractionExpansion
  const site_t n_site;                                //number of sites
  const frequency_t n_matsubara;        //number of matsubara freq
  const frequency_t n_matsubara_measurements;        //number of measured matsubara freq
  const itime_index_t n_tau;                        //number of imag time slices
  const itime_t n_tau_inv;                        //the inverse of n_tau
  const frequency_t n_self;                        //number of self energy (W) binning points
  const boost::uint64_t mc_steps;
  const std::uint64_t therm_steps;

  const double beta;
  const double temperature;                        //only for performance reasons: avoid 1/beta computations where possible
  const double alpha;
  const U_matrix U;
  std::vector<std::pair<spin_t, spin_t>> interaction_pairs;


  const unsigned int recalc_period;
  const unsigned int measurement_period;

  /*InteractionExpansion's roundoff threshold*/
  const double almost_zero;

  /*private member variables*/
  itime_green_function_t bare_green_itime;

  vertex_array vertices;
  std::vector<inverse_m_matrix> M;

  double sign;
  unsigned int measurement_method;


  std::uint64_t step;

};



#endif
