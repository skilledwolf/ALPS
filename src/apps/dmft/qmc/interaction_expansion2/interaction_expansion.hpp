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

#include <alps/ngs.hpp>
#include <alps/mcbase.hpp>

#include <alps/alea.h>
#include <cmath>
#include "green_function.h"
#include "alps_solver.h"
#include "types.h"
#include "solver.h"
#include "alps_solver.h"
#include "fouriertransform.h"
#include "U_matrix.h"
#include "../ctint/operator.hpp"
#include "../ctint/green_matrix.hpp"
#include <alps/numeric/matrix.hpp>


#include "../ctint/types.hpp"

namespace alps { namespace ctint { struct update_kernel; } }


/*class InteractionExpansionSim: public alps::scheduler::MCSimulation, public alps::MatsubaraImpurityTask
{
public:

  InteractionExpansionSim(const alps::ProcessList &w, const boost::filesystem::path &p) : alps::scheduler::MCSimulation(w,p) {}
  
  InteractionExpansionSim(const alps::ProcessList &w, const alps::Parameters &p) : alps::scheduler::MCSimulation(w,p) {p_=p;}
  
  std::pair<matsubara_green_function_t,itime_green_function_t> get_result(); 

  void evaluate_selfenergy_measurement_matsubara(const alps::ObservableSet &gathered_measurements, 
                                                 matsubara_green_function_t &green_matsubara_measured,
                                                 const matsubara_green_function_t &bare_green_matsubara, 
                                                 std::vector<double>& densities, const double &beta, 
                                                 const int n_site, const int n_flavors, const int n_matsubara) const;

  void evaluate_selfenergy_measurement_itime_rs(const alps::ObservableSet &gathered_measurements, itime_green_function_t &green_result,
                                                const itime_green_function_t &green0, const double &beta, const int n_site, 
                                                const int n_flavors, const int n_tau, const int n_self) const;

  double green0_spline(const itime_green_function_t &green0, const itime_t delta_t, const int s1, const int s2, 
                       const spin_t flavor, int n_tau, double beta) const;
  
private:

  alps::Parameters p_;
};*/



class InteractionExpansion: public alps::mcbase
{
  friend struct alps::ctint::update_kernel;
  double ctint_uniform() { return random(); }
  unsigned int ctint_site(unsigned int sites) { return static_cast<unsigned int>(random() * sites); }
  [[noreturn]] void ctint_invalid_flavors(const char* message) { throw std::logic_error(message); }

public:

  InteractionExpansion(const alps::params& p, int rank);
  ~InteractionExpansion() {}
  bool is_thermalized() const {return true;} //thermalization is done in the constructor. It's not a big deal here.
  void update();
  void measure();
  double fraction_completed() const;
    
protected:
  
  /*functions*/
  /*io & initialization*/
  void initialize_simulation(const alps::params &parms); // called by constructor
  // in file io.cpp
  void read_bare_green(std::ifstream &G0_omega, std::ifstream &G0_tau);
  void print(std::ostream &os);
  
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
  
  /*abstract virtual functions. Implement these for specific models.*/
  virtual double try_add()=0;
  virtual void perform_add()=0;
  virtual void reject_add()=0;
  virtual double try_remove(unsigned int vertex_nr)=0;
  virtual void perform_remove(unsigned int vertex_nr)=0;
  virtual void reject_remove()=0;
  
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
  const unsigned long therm_steps;                
  const double max_time_in_seconds;
  
  const double beta;                                
  const double temperature;                        //only for performance reasons: avoid 1/beta computations where possible        
  const double onsite_U;                        
  const double alpha;                                
  const U_matrix U;
  
  
  const unsigned int recalc_period;                
  const unsigned int measurement_period;        
  const unsigned int convergence_check_period;        
  
  /*InteractionExpansion's roundoff threshold*/
  const double almost_zero;                        
  /*PRNG seed*/
  const int seed;                                
  
  /*private member variables*/
  matsubara_green_function_t green_matsubara;
  matsubara_green_function_t bare_green_matsubara;
  itime_green_function_t bare_green_itime;
  itime_green_function_t green_itime;
  std::vector<green_matrix> g0;
  boost::shared_ptr<FourierTransformer> fourier_ptr;
  
  vertex_array vertices;
  std::vector<inverse_m_matrix> M;
    
  double weight;
  double sign;
  unsigned int measurement_method;
  bool thermalized;
  
  simple_hist pert_hist;
  unsigned int hist_max_index;
  simple_hist **vertex_histograms;
  unsigned int vertex_histogram_size;
  
  unsigned long step;        
  time_t start_time;
  clock_t update_time;
  clock_t measurement_time;

};



/*aux functions*/
std::ostream& operator << (std::ostream& os, const std::vector<double>& v);
std::ostream& operator << (std::ostream &os, const vertex_array &vertices);
std::ostream& operator << (std::ostream &os, const vertex &v);
std::ostream& operator << (std::ostream &os, const c_or_cdagger &c);
std::ostream& operator << (std::ostream& os, const simple_hist &h);



//Use this for the most simple single site Hubbard model.
class HalfFillingHubbardInteractionExpansion: public InteractionExpansion{
public:
  HalfFillingHubbardInteractionExpansion(const alps::params& p, int rank)
    :InteractionExpansion(p, rank)
  {
    if(n_flavors !=1){throw std::invalid_argument("you need a different model for n_flavors!=1.");}
  }
  double try_add();
  void perform_add();
  void reject_add();
  double try_remove(unsigned int vertex_nr);
  void perform_remove(unsigned int vertex_nr);
  void reject_remove();
};



class HubbardInteractionExpansion: public InteractionExpansion{
public:
  HubbardInteractionExpansion(const alps::params& p, int rank)
    :InteractionExpansion(p, rank)
  {
    if(n_flavors !=2){throw std::invalid_argument("you need a different model for n_flavors!=2.");}
  }
  double try_add();
  void perform_add();
  void reject_add();
  double try_remove(unsigned int vertex_nr);
  void perform_remove(unsigned int vertex_nr);
  void reject_remove();
};



//Use this for multiple bands where you have terms Un_i n_j
class MultiBandDensityHubbardInteractionExpansion: public InteractionExpansion{
public:
  MultiBandDensityHubbardInteractionExpansion(const alps::params& p, int rank)
    :InteractionExpansion(p, rank)
  {
    if(n_site !=1){throw std::invalid_argument("you need a different model for n_site!=1.");}
  }
  
  double try_add();
  void perform_add();
  void reject_add();
  double try_remove(unsigned int vertex_nr);
  void perform_remove(unsigned int vertex_nr);
  void reject_remove();
};

#endif
