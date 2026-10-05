/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2001-2009 by Matthias Troyer <troyer@comp-phys.org>,
*                            Simon Trebst <trebst@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#ifndef ALPS_APPLICATIONS_WORM_RUN_H
#define ALPS_APPLICATIONS_WORM_RUN_H

#define P_REMOVE 0.5


#include <alps/lattice.h>
#include <alps/expression.h>

#include "WKink.h"
#include "random.h"
#include "../qmc.h"


#include <algorithm>
#include <alps/hdf5/stdarray.hpp>

#include <boost/multi_array.hpp>

//#define Print_config
//#define Print_steps
//#define Print_spins
//#define Print_update
//#define Print_detail
//#define Print_weight

//#define CHECK_OFTEN 

//#define SIMPLE

// non-local interactions?

// data structure for kinks, default is list
// #define USE_VECTOR
// #define USE_SET

struct subinterval_info{ time_struct start_time;
                         time_struct end_time;
                         double delta_t;
                         double delta_e;
                         double integrated_time;
                         double integrated_weight; };

struct worm_state {
  using state_type=uint8_t;
  using kink_type=Kink<state_type>;
  using kinklist_type=std::list<kink_type>;
  using kinkvector_type=std::vector<kinklist_type>;
  uint64_t steps=0;
  int measurements_done=1, num_kinks=0;
  double worms_per_update=1., Sign=1.;
  unsigned last_id_=0;
  int corrections_upwards=0, corrections_downwards=0;
  bool preadjustment_done=false, adjustment_done=false;
  std::vector<int> nob;
  kinkvector_type kinks;
  std::vector<state_type> initial_state_;
};

class WRun : public QMCRun<>, private worm_state {
public:

  static void print_copyright(std::ostream&);
  WRun(alps::params const&,size_t bins=128,size_t chain=0);
  alps::params sampling_parameters() const {
    auto result=parameters;
    if (canonical) result[adjust_parameter]=double(parms[adjust_parameter]);
    return result;
  }
  void save(alps::hdf5::archive&) const override;
  void load(alps::hdf5::archive&) override;
  void update() override { dostep(); }
  void measure() override {} // Measurement interval handled by dostep.
  double fraction_completed() const override { return work_done(); }
  uint64_t completed_sweeps() const { return steps; }

  void dostep();
  bool is_thermalized() const;
  double work_done() const;

private:

  //- types --------------------------------------------------

  // model specification
  int number_of_bond_types; 
  std::vector<uint8_t> site_state; 

  std::vector<std::pair<uint8_t,uint8_t> > site_type_for_bond_type_;
  std::vector<int> original_bond_type;
  
  std::vector< std::vector<double> > matrix_element_raise_;
  std::vector< std::vector<double> > matrix_element_lower_;

  std::vector< std::vector<double> > site_matrix, requested_site_matrix;
  double adjustment_energy_shift=0.;
  std::vector<double>                hopping_matrix;

  std::map<int, boost::multi_array<double,2> > diagonal_matrix_element;
  
  using state_type=worm_state::state_type;
  size_t chain_;
  typedef ::cyclic_iterator<kinklist_type> cyclic_iterator;
  typedef cyclic_iterator::base_iterator iterator;

  typedef Wormhead<kinklist_type,graph_type> wormhead_type;

  //- member functions ---------------------------------------

  // worm creation/annihilation
  inline bool create_worm();
  inline bool annihilate_worm();
  int64_t make_worm();

  // update methods
  void shift_kink(wormhead_type&,wormhead_type&);
  void insert_jump(wormhead_type&,wormhead_type&,int,int);
  void remove_jump(wormhead_type&,wormhead_type&,int,int);

  // unperturbed energy of site states, includes neighbor terms
  inline double H0(const site_descriptor&);

  // unperturbed energy differences
  inline double Delta_H0(const site_descriptor&, const state_type&, const state_type&, std::vector<state_type>&);
  inline double Delta_H0(const site_descriptor&, const state_type&, const state_type&, std::vector<state_type>&,
                         const site_descriptor&, const state_type&, const state_type&, std::vector<state_type>&, const bond_descriptor&);

  inline void Update_Delta_H0(double& delta_e, const state_type& state1, const state_type& state2,  
                              const state_type& nb_state1, const state_type& nb_state2, 
                              const bond_descriptor& b) {
    if (nonlocal) {
      delta_e -= neighbor_energy(state1, nb_state1, b);
      delta_e += neighbor_energy(state2, nb_state1, b);
      delta_e += neighbor_energy(state1, nb_state2, b);
      delta_e -= neighbor_energy(state2, nb_state2, b);
    }
  }   // Update_Delta_H0

  // calculate integrated weight taking into account all subintervals
  inline double integrated_weight(const double& lambda, const double& time)
    {
#ifdef SIMPLE
      if(lambda!=0.)
        return 0.;
      else
        return time;
#else
      if(lambda==0. || fabs(lambda*time)<1e-10)
        return time;
      if ((-lambda*time > log_numeric_limits_double ) && is_thermalized())
        boost::throw_exception(std::logic_error("Exceeding std::numeric_limits<double>::max() in integrated weight"));        
      return 1./lambda*(1.-std::exp(-lambda*time));
#endif
    }   // WRun::integrated_weight

  // explore all subintervals due to non-local interaction
  void traverse_subintervals(wormhead_type&);
  void traverse_subintervals(const time_struct&, const time_struct&, const double&, 
                             const site_descriptor&,
                             const state_type&, const state_type&, double);

  // find time steps where neighbor states change
  std::pair<time_struct, time_struct> adjacent_subinterval(const site_descriptor&, const time_struct&,
                                                                      std::vector<state_type>&);

  // calculate worm creation probability taking into account all subintervals
  double worm_creation_probability(const time_struct&, const double&, const site_descriptor&,
                                   const state_type&, const state_type&, double);
                                    
  // measurements
  void make_meas();
  inline void sample() {
    if(--measurements_done==0) {
      measurements_done=skip_measurements;
      make_meas();
    }
  }

  // adjustment
  int get_particle_number();
  void adjustment();
  void set_adjusted_parameter(double);
  bool canonical;
  std::string adjust_parameter;

  // checks
  void print_spins();
  void check_spins();

  void start();
  
  //- model related functions --------------------------------

  state_type min_state() { return std::numeric_limits<state_type>::min(); }
  state_type max_state() { return std::numeric_limits<state_type>::max(); }

  // state_type initial_state()                 { return state_type(0); }
  state_type create(state_type s, int=0)     { return s+1; }
  state_type annihilate(state_type s, int=0) { return s-1; }

  inline double creation_matrix_element(const state_type& state, const bool& c, site_descriptor site) { 
    int this_site_type = site_type(site);
    return (c ? matrix_element_raise_[this_site_type][state] * matrix_element_raise_[this_site_type][state] 
            : matrix_element_lower_[this_site_type][state] * matrix_element_lower_[this_site_type][state]);
  }

  inline double onsite_energy(const state_type& state, site_descriptor site) { 
    return site_matrix[inhomogeneous_site_type(site)][state];
  }

  inline double neighbor_energy(const state_type& state1, const state_type& state2, const bond_descriptor& b) {
    return diagonal_matrix_element[bond_type[b]][state1][state2];
  }   

  inline double hopping_matrix_element(const bond_descriptor& b) {
    return hopping_matrix[bond_type[b]];
  }

  // initialization
  void initialize_hamiltonian();
  std::vector<double> hamiltonian_state(bool include_sites=true) const;
  boost::multi_array<double,4> bond_hamiltonian(const bond_descriptor&);
  std::vector<double> site_hamiltonian(const site_descriptor&);
  void print_hamiltonian();
  void create_observables();

  //- data ---------------------------------------------------

  std::vector<wormhead_type> worm_head;
  std::valarray<double> stat;
  double eta;
  uint64_t thermal_sweeps;
  int skip_measurements;
  bool have_worm;
  bool chain_kappa;
  int num_chains;

  int worms_per_kink;
  double log_numeric_limits_double;
  
  std::vector<subinterval_info> subinterval;
  bool subinterval_valid;
  int current_head_num;

  bool nonlocal;
  
  bool use_1D_stiffness ; //@#$br


  //- kinks and related data structures ----------------------

  // iterator to first kink (in time)
  cyclic_iterator first_kink(site_descriptor i) { return cyclic_iterator(kinks[i],kinks[i].begin());} 
  state_type initial_state(site_descriptor i) 
  { return kinks[i].empty() ?  initial_state_[i] : kinks[i].rbegin()->state(); }
  std::vector<int> chain_number, chain_size;
    
  alps::property_map<alps::bond_type_t,graph_type,int>::type bond_type;

  inline void erase_kink(int site, iterator w) 
  {
    kinks[site].erase(w);
  }

  inline iterator insert_kink(int site, iterator w, kink_type kink)
  {
    kinklist_type& l(kinks[site]);
#ifdef USE_SET
      w=l.insert(w,kink);
#else
    if (l.empty() || kink.time() < l.begin()->time())
      w=l.insert(l.begin(),kink);
    else if (w==l.begin())
      w=l.insert(l.end(),kink);
    else
      w=l.insert(w,kink);
#endif
    return w;
  }

  inline iterator move_kink(int site, iterator w, time_struct newtime)
  {
    kinklist_type& l(kinks[site]);
#ifdef USE_SET
    kink_type kink(*w);
    l.erase(w);
    kink.set_time(newtime);
    w=l.insert(kink).first;
#else
    iterator k=w;
    ++k;
    if (w==l.begin() && newtime > l.rbegin()->time()) {
      kink_type kink(*w);
      l.erase(w);
      w=l.insert(l.end(),kink);
    }
    else if (k==l.end() && newtime < l.begin()->time()) {
      kink_type kink(*w);
      l.erase(w);
      w=l.insert(l.begin(),kink);
    }
    w->set_time(newtime);
#endif
    return w;
  }
};

//- Unperturbed energy ----------------------------------------------------------

inline double WRun::H0(const site_descriptor& s1) {
  //
  //  determines the energy of the state on site s1 at time 0 taking into account 
  //  the states of all neighbors at that particular time.
  //  This function does not doublecount the bond energy terms.
  //

  // determine onsite energy
  state_type state1 = initial_state(s1);
  double Result = onsite_energy(state1,s1);

  if (nonlocal) {
    // determine neighbor energy
    state_type state2;
    neighbor_bond_iterator nbi, nbi_end;
    int nb=0;
    for(boost::tie(nbi, nbi_end) = neighbor_bonds(s1); nbi != nbi_end; ++nbi) {
      state2 = initial_state(neighbor(s1, nb));
      Result += 0.5*neighbor_energy(state1, state2,*nbi);
      nb++;
    }  
  }

  return Result;
}   // WRun::H0


#endif
