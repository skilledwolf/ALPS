/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1999-2009 by Matthias Troyer <troyer@comp-phys.org>,
*                            Fabian Stoeckli <fabstoec@student.ethz.ch>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#ifndef ALPS_APPLICATIONS_MC_SPIN_ABSTRACTSPINMC_H_
#define ALPS_APPLICATIONS_MC_SPIN_ABSTRACTSPINMC_H_

#include <alps/scheduler.h>
#include <alps/alea.h>
#include <boost/property_map/vector_property_map.hpp>

#include <alps/model/sign.h>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include "couplings.h"
#include "matrices.h"
#include "settings.h"

typedef alps::uint64_t uint64_t;

struct update_info_type 
{
  std::size_t clustersize;
  double m2_;
};

 
template<class MAT>
class AbstractSpinSim : public alps::scheduler::LatticeMCRun<> {
public : 
  AbstractSpinSim(const alps::ProcessList&,const alps::Parameters&,int);
  ~AbstractSpinSim();  
 
  bool change_parameter(const std::string& name, const alps::StringValue& value);
  void dostep();
  void save(alps::ODump& dump) const;
  void load(alps::IDump& dump);
  bool is_thermalized() const;
  double work_done() const;

  virtual update_info_type do_update()=0;
  virtual void do_measurements(update_info_type)=0;

protected: 
  // Resolved first, directly after the base class has built the lattice.
  spinmc::settings<MAT::dim> settings_;
  double beta_;
  uint64_t sweeps_;
  TinyVector<double,MAT::dim> h_;
  TinyVector<double,MAT::dim> h_normalized;
  bool has_magnetic_field_;
  uint64_t sweeps_done_;
  uint64_t thermalization_;
  double thermalization_fraction_;
  uint64_t thermalization_sweeps_;
  
  bool cluster_updates_;
  bool general_case_;
  typedef boost::property_map<graph_type,alps::bond_type_t>::type  bond_type_map;
  typedef boost::property_map<graph_type,alps::site_type_t>::type  site_type_map;
  alps::property_map<alps::bond_type_t,graph_type,int>::type bond_type;

  boost::vector_property_map<double,site_type_map> spinfactor_; /* S */
  boost::vector_property_map<MAT,site_type_map> selfinteraction_; /* D */
  boost::vector_property_map<MAT,bond_type_map> couplings_; /* J */ 
  boost::vector_property_map<double,bond_type_map> couplings_det; /* |J| */
 
  bool use_error_limit;
  bool ferromagnetic_;
  bool antiferromagnetic_;
  bool quantum_convention_;
  double g_;
  int print_sweeps_;
};  

template<class MAT>
AbstractSpinSim<MAT>::AbstractSpinSim(const alps::ProcessList& w, 
        const alps::Parameters& myparms,int n)
  : alps::scheduler::LatticeMCRun<>(w,myparms,n),
  settings_(spinmc::resolve_settings<MAT::dim>(parms)),
  beta_(settings_.beta),
  sweeps_(settings_.sweeps),
  h_(settings_.h),
  sweeps_done_(0),
  thermalization_(settings_.thermalization),
  thermalization_fraction_(0.),
  thermalization_sweeps_(0),
  general_case_(false),
  bond_type(alps::get_or_default(alps::bond_type_t(),graph(),0)),
  spinfactor_(boost::get(alps::site_type_t(),graph())),
  selfinteraction_(boost::get(alps::site_type_t(),graph())),
  couplings_(boost::get(alps::bond_type_t(),graph())),
  couplings_det(boost::get(alps::bond_type_t(),graph())),
  ferromagnetic_(true),
  antiferromagnetic_(true),
  quantum_convention_(settings_.quantum_convention),
  g_(settings_.g),
  print_sweeps_(settings_.print_sweeps)
{
  // 2009-01-27: wistaria@comp-phys.org
  // moved from the inilializaton list to avoid 'internal compiler error' in gcc 4.3.x
  h_normalized = h_*(1./std::sqrt(h_.get_length2()));
  has_magnetic_field_ = h_normalized.get_length2()>1.e-10;

  if (inhomogeneous())
    boost::throw_exception(std::runtime_error("Disordered lattices not supported by the classical Monte Carlo program.\n"));

  use_error_limit = settings_.use_error_limit;

  if (has_magnetic_field_ && quantum_convention_) 
    boost::throw_exception(
      std::runtime_error("Cannot use quantum convention with nonzero magnetic field ."));


  // Reading UPDATE through the worker's mutable parameters has always
  // inserted an empty value when it is absent, and that entry is saved with
  // the run. Keep the insertion until the saved parameter set may change.
  (void)parms["UPDATE"];

  if (settings_.update==spinmc::update_request::cluster && has_magnetic_field_)
    {boost::throw_exception(std::runtime_error("Illegal update type " 
               + settings_.update_text
               + " with magnetic field not set to zero."));}

  // Site types in the order of the sites, and bond classes (bond type,
  // source and target site type) in the order of the bonds. Couplings are
  // resolved per class, and the class index becomes the bond type.
  std::vector<int> site_types;
  std::map <int,bool> site_type_visited;
  for (site_iterator si = sites().first; si!=sites().second; si++) {
    int type = site_type(*si);
    if (!site_type_visited[type]) {
      site_type_visited[type]=true;
      site_types.push_back(type);
    }
  }

  typedef boost::tuple<int,int,int> bond_tuple_type;
  std::map<bond_tuple_type,int> new_bond_type;
  std::vector<spinmc::bond_class> bond_classes;
  for (bond_iterator b = bonds().first; b!=bonds().second ;++b) {
    bond_tuple_type this_bond(bond_type[*b],site_type(source(*b)),site_type(target(*b)));
    if (new_bond_type.find(this_bond)==new_bond_type.end()) {
      new_bond_type[this_bond] = bond_classes.size();
      spinmc::bond_class c = {this_bond.get<0>(),this_bond.get<1>(),this_bond.get<2>()};
      bond_classes.push_back(c);
    }
  }

  const spinmc::couplings<MAT> terms = spinmc::resolve_couplings<MAT>(
      parms, site_types, bond_classes, quantum_convention_, settings_.default_coupling);
  general_case_ = terms.general_case;
  ferromagnetic_ = terms.ferromagnetic;
  antiferromagnetic_ = terms.antiferromagnetic;

  std::size_t next_site_type = 0;
  site_type_visited.clear();
  for (site_iterator si = sites().first; si!=sites().second; si++) {
    int type = site_type(*si);
    if (!site_type_visited[type]) {
      site_type_visited[type]=true;
      spinfactor_[*si]=terms.sites[next_site_type].spinfactor;
      selfinteraction_[*si]=terms.sites[next_site_type].selfinteraction;
      ++next_site_type;
    }
  }

  std::vector<bool> coupling_set(bond_classes.size(), false);
  for (bond_iterator b = bonds().first; b!=bonds().second ;++b) {
    bond_tuple_type this_bond(bond_type[*b],site_type(source(*b)),site_type(target(*b)));
    int num_bond_type = new_bond_type[this_bond];
    bond_type[*b] = num_bond_type;
    if (!coupling_set[num_bond_type]) {
      couplings_[*b]=terms.bonds[num_bond_type].coupling;
      couplings_det[*b] = terms.bonds[num_bond_type].det;
      coupling_set[num_bond_type] = true;
    }
  }

  if (settings_.update==spinmc::update_request::cluster) {
    if (MAT::allow_cluster_update)
      cluster_updates_=true;
    else {
      std::cerr << "Invalid input: \n"
                << "Cluster updates are not allowed with the given matrix "
                << "types.\nCluster updates are only allowed for matrices "
                << "reducing to a\nsingle element (these are diagonal matrices "
                << "with always the same\nelement on the diagonal). For this "
                << "case, give one single value\nin the parameter file.\n";
      boost::throw_exception(std::runtime_error("wrong parameter submitted"));
    }
  }
  else if (settings_.update==spinmc::update_request::local)
    cluster_updates_=false;
  else if (settings_.update==spinmc::update_request::automatic)
    cluster_updates_ = !alps::is_frustrated(graph(),couplings_det);
  else
    boost::throw_exception(std::runtime_error("Illegal update type " + settings_.update_text));

  cluster_updates_ = cluster_updates_ && MAT::allow_cluster_update;

  if ( cluster_updates_) std::cerr << "Using cluster updates\n";
  // check whether frustrated 
  measurements << alps::RealObservable("Energy");
  measurements << alps::RealVectorObservable("Bond-type Energy");
  measurements << alps::RealObservable("Energy Density");
  measurements << alps::RealVectorObservable("Bond-type Energy Density");
  measurements << alps::RealObservable("Energy^2");
  measurements << alps::RealObservable("beta * Energy / sqrt(N)");
  measurements << alps::RealObservable("(beta * Energy)^2 / N");
  measurements << alps::RealObservable("|Magnetization|");
  measurements << alps::RealObservable("Magnetization along Field");
  measurements << alps::RealObservable("Magnetization^2");
  measurements << alps::RealObservable("E.Magnetization^2");
  measurements << alps::RealObservable("Magnetization^4");
  measurements << alps::RealObservable("E.Magnetization^4");
  measurements << alps::RealObservable("Susceptibility");
  if(cluster_updates_)
    measurements << alps::RealObservable("Cluster size");
  if (is_bipartite()) {
    measurements << alps::RealObservable("|Staggered Magnetization|");
    measurements << alps::RealObservable("Staggered Magnetization^2");
  }
  
}



template<class MAT>
AbstractSpinSim<MAT>::~AbstractSpinSim() 
{
// nothing to be done
}  
 
// Check if simulation is finished
// take the measurement that is closer to be finished
// ie. the amount of work done is max(work_done(sweeps),work_done(error_limit)).
template<class MAT>
double AbstractSpinSim<MAT>::work_done() const {
  if (!is_thermalized())
    return 0.;
  double wd_error = 0.0;
  double wd_sweeps = 0.0;
  if (use_error_limit) {
    alps::RealObservable myObs;
#ifndef BOOST_NO_EXCEPTIONS
    try {
#endif
      myObs = measurements.template get<alps::RealObservable>(parms["ERROR_VARIABLE"].c_str());
#ifndef BOOST_NO_EXCEPTIONS
    }
    catch (std::exception& e) {
      std::cerr << "Caught exception: " << e.what() << "\n";
      std::exit(-5);
   }
#endif
    if (myObs.count() == 0)
      return 0.;

    double current_error = myObs.error();
    double error_limit = parms["ERROR_LIMIT"];
    double tmp = error_limit/current_error;
    wd_error = tmp*tmp;
  }
  double good_sweeps = sweeps_done_-thermalization_sweeps_;
  wd_sweeps = (good_sweeps)/double(sweeps_);

  if (wd_sweeps >= 1.0)
    return wd_sweeps;

  if (wd_error > wd_sweeps) return wd_error; 
  else return wd_sweeps;
}

// It is desirable to be able to increase the number of SWEEPS
template<class MAT>
bool AbstractSpinSim<MAT>::change_parameter(const std::string& name, const alps::StringValue& value) {
  uint64_t new_sweeps = 0;
  if(name=="SWEEPS")
    new_sweeps = uint64_t(value);
  // Is it sensible to do it ?
  if(new_sweeps > 0) {
    sweeps_ = new_sweeps;
    return true;
   }
  // Otherwise we cannot do it
  return false;
}

// Do a MC step - and do measurements (if thermalized)
// Do a MC step only sweeps_ times.
template<class MAT>
void AbstractSpinSim<MAT>::dostep() {
  ++sweeps_done_;
  update_info_type s= do_update();
  if (!is_thermalized()) {
    ++thermalization_sweeps_;
    thermalization_fraction_+=double(s.clustersize)/double(num_sites());
  }
  else {
// Is there any good reason to discard sweeps beyond the required ones ?
// if (sweeps_done_ < sweeps_) {
    do_measurements(s);
  }
}
  
// Check if simulation is thermalized
template<class MAT>
bool AbstractSpinSim<MAT>::is_thermalized() const {
  return (thermalization_fraction_ >= thermalization_); 
}


// save simulation
template<class MAT>
void AbstractSpinSim<MAT>::save(alps::ODump& dump) const {
  dump << sweeps_done_ << thermalization_fraction_ << thermalization_sweeps_;
}

// load simulation
template<class MAT>
void AbstractSpinSim<MAT>::load(alps::IDump& dump) {
  dump >> sweeps_done_ >> thermalization_fraction_ >> thermalization_sweeps_;
  if(where.empty())
    measurements.compact();
}

#endif
