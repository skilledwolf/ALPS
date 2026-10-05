/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2001-2006 by Fabien Alet <alet@comp-phys.org>,
*                            Matthias Troyer <troyer@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#ifndef ALPS_APPLICATIONS_QMC_H
#define ALPS_APPLICATIONS_QMC_H

#include <boost/math/constants/constants.hpp>
#include <alps/mcbase.hpp>
#include <alps/model.h>
#include <alps/ngs/make_deprecated_parameters.hpp>
#include "simulation.hpp"
#include <alps/scheduler/measurement_operators.h>
#include <set>

#include <boost/optional.hpp>
#include <boost/assert.hpp>

template <class G=typename alps::graph_helper<>::graph_type, class StateType = boost::uint8_t>
class QMCRun : public native_qmc::simulation, public alps::graph_helper<G>, public alps::model_helper<>
             , public alps::MeasurementOperators
{
public :
  typedef StateType state_type;
  typedef alps::graph_helper<G> super_type;
  typedef typename super_type::site_iterator site_iterator;

  QMCRun(alps::params const&,size_t,size_t,bool=false);
  unsigned winding_dimension() const { return winding_dimension_; }
  double density_reference() const { return density_reference_; }
  void save(alps::hdf5::archive& ar) const override {
    alps::mcbase::save(ar);
    ar["checkpoint/density_reference"] << density_reference_;
  }
  void load(alps::hdf5::archive& ar) override {
    double reference;
    ar["checkpoint/density_reference"] >> reference;
    alps::mc::batch density;
    if (is_charge_model_) {
      alps::alea::hdf5_serializer codec(ar,"measurements");
      alps::alea::deserialize(codec,"Centered Density Moments",density);
    }
    if (std::isinf(reference) || (is_charge_model_ && density.count() && !std::isfinite(reference)))
      throw std::invalid_argument("Invalid density reference");
    alps::mcbase::load(ar);
    density_reference_=reference;
  }
  template<class T> void record(std::string const& name,T const& value,double sign) {
    if (!recording_) return;
    if (name=="Density") {
      auto adapter=alps::alea::make_adapter(value);
      alps::alea::column<double> sample=alps::alea::column<double>::Zero(adapter.size());
      adapter.add_to(alps::alea::view<double>(sample.data(),adapter.size()));
      double density=sample[0]*sign;
      if (std::isnan(density_reference_)) density_reference_=density;
      double delta=density-density_reference_;
      alps::mc::record(*this,"Centered Density Moments",alps::alea::column<double>{delta*sign,delta*delta*sign,sign});
    }
    native_qmc::simulation::record(name,value,sign);
  }
  void record(std::string const& name,std::valarray<double> const& value,double sign) {
    record(name,std::vector<double>(std::begin(value),std::end(value)),sign);
  }
  int random_int(int n) { return int(random()*n); }
  double random_real() const { return random(); }
  double random_01() const { return random(); }
  bool has_sign_problem() const {
    return alps::has_sign_problem(this->model(),static_cast<super_type const&>(*this),parms);
  }
protected:
  static alps::Parameters graph_parameters(alps::params const& p) {
    alps::Disorder::seed(p.value_or<uint32_t>("DISORDER_SEED",uint32_t(p.value_or("SEED",42))));
    return alps::make_deprecated_parameters(p);
  }
  unsigned winding_dimension_=0;
  alps::Parameters parms;
  double density_reference_=std::numeric_limits<double>::quiet_NaN();
protected:
  double beta;
  bool is_spin_model_;
  bool is_charge_model_;
  bool measure_local_density_;
  bool measure_local_magnetization_;
  bool measure_site_type_density_;
  bool measure_local_compressibility_;
  bool measure_site_compressibility_;
  bool measure_correlations_;
  bool measure_structure_factor_;
  bool measure_green_function_;
  bool measure_bond_type_stiffness_;
  bool measure_CBS_order_;
  boost::optional<typename super_type::site_descriptor> measurement_origin_;
  std::vector<unsigned int> distance_mult;
  std::size_t num_site_types_;
  std::size_t num_bond_types_;
  std::valarray<double> site_type_density_;
  std::valarray<double> local;

  unsigned maximum_sitetype;
  state_type maximum_number_of_states;
  std::map<int,int> number_states_for_site_type_;
  std::vector<int>        site_site_type;
  std::vector<state_type> site_number_of_states;

  std::map<std::string,std::vector<std::vector<double> > > diagonal_matrix_element;

  boost::optional<int> restricted_particle_number;
  boost::optional<double> restricted_magnetization;


  void create_common_observables();
  bool do_common_measurements(double sign, const std::vector<state_type>& local,
            const std::valarray<double>& local_int=std::valarray<double>());
  void initialize_site_states();
private:
  bool build_diagonal_operator(std::string const& op);
};

template <class G, class StateType>
QMCRun<G,StateType>::QMCRun(alps::params const& p,size_t bins,size_t chain,bool issymbolic)
  : native_qmc::simulation(p,bins,chain), super_type(graph_parameters(p)),
    alps::model_helper<>(static_cast<super_type const&>(*this),alps::make_deprecated_parameters(p),issymbolic),
    alps::MeasurementOperators(alps::make_deprecated_parameters(p)),
    parms([&] {
      auto values=alps::make_deprecated_parameters(p);
      values.copy_undefined(this->model().default_parameters());
      return values;
    }())
  , beta(parms.defined("Beta") ? alps::evaluate<double>(parms["Beta"],parms)
      :  (parms.defined("beta") ? alps::evaluate<double>(parms["beta"],parms)
        : (parms.defined("BETA") ? alps::evaluate<double>(parms["BETA"],parms)
           : (parms.defined("T") ? 1./alps::evaluate<double>(parms["T"],parms)
              : (parms.defined("TEMPERATURE") ? 1./alps::evaluate<double>(parms["TEMPERATURE"],parms)
                 : 1./alps::evaluate<double>(parms["temperature"],parms)))))),
    is_spin_model_(false),
    is_charge_model_(false),
    measure_local_density_(false),
    measure_local_magnetization_(false),
    measure_site_type_density_(parms.value_or_default("MEASURE[Site Type Density]",false)),
    measure_local_compressibility_(false),
    measure_site_compressibility_(false),
    measure_correlations_(parms.value_or_default("MEASURE[Correlations]",false)),
    measure_structure_factor_(parms.value_or_default("MEASURE[Structure Factor]",false)),
    measure_green_function_(parms.value_or_default("MEASURE[Green Function]",false)),
    measure_bond_type_stiffness_(parms.value_or_default("MEASURE[Bond Type Stiffness]",false)),
    measure_CBS_order_(parms.value_or_default("MEASURE[CBS Order]",false)),
    num_site_types_(alps::maximum_vertex_type(this->graph())+1),
    num_bond_types_(alps::maximum_edge_type(this->graph())+1)

{
  is_signed_=has_sign_problem();
  if (!this->num_sites() || !this->num_bonds() || bins<2 || bins%2)
    throw std::invalid_argument("QMC requires a nonempty lattice and even batch capacity >= 2");
  for (auto [it,end]=this->bonds();it!=end;++it)
    if (this->source(*it)==this->target(*it)) throw std::invalid_argument("QMC does not support self-bonds");
  if (parms.defined("INITIAL_SITE"))
    measurement_origin_=static_cast<int>(parms["INITIAL_SITE"]);
  if (measurement_origin_ && *measurement_origin_>=this->num_sites())
    throw std::invalid_argument("INITIAL_SITE is outside the lattice");
  winding_dimension_=this->dimension();
  for (auto [it,end]=this->bonds();it!=end;++it) {
    auto const& v=this->bond_vector_relative(*it);
    if (v.size()!=winding_dimension_ || !std::all_of(v.begin(),v.end(),[](double x){return std::isfinite(x);}))
      winding_dimension_=0;
  }
  if (measure_CBS_order_ && this->dimension()<2)
    throw std::invalid_argument("CBS Order requires at least two dimensions");
  if (measure_bond_type_stiffness_)
    throw std::invalid_argument("Bond-type stiffness is not implemented by these QMC kernels");
  if (parms.defined("RESTRICT_MEASUREMENTS[N]"))
    restricted_particle_number.reset(static_cast<int>(parms["RESTRICT_MEASUREMENTS[N]"]));
  if (parms.defined("RESTRICT_MEASUREMENTS[Sz]"))
    restricted_magnetization.reset(static_cast<double>(parms["RESTRICT_MEASUREMENTS[Sz]"]));


  if(!std::isfinite(beta) || beta<=0)
    boost::throw_exception(
         std::out_of_range("inverse temperature beta must be positive and finite"));

}

template <class G, class StateType>
bool QMCRun<G,StateType>::build_diagonal_operator(std::string const& name)
{
  if (diagonal_matrix_element.find(name) != diagonal_matrix_element.end())
    return true;
  std::vector<std::vector<double> > vec(maximum_sitetype + 1);

  std::map<int,int>::const_iterator it = number_states_for_site_type_.begin();
  for (; it != number_states_for_site_type_.end(); ++it) {
    int sitetype = it->first;
    int number_states = it->second;
    if (number_states) {
      alps::SiteOperator term(name);
      boost::multi_array<alps::Expression,2> matrix_symbolic =
          alps::get_matrix(alps::Expression(),term,this->model().basis().site_basis(sitetype),this->parms);
      for (int i=0;i<number_states;++i)
        for (int j=0;j<number_states;++j) {
          if (!matrix_symbolic[i][j].can_evaluate())
            return false;
          else if (i!=j && alps::evaluate<double>(matrix_symbolic[i][j]))
            return false;
          else if (i==j)
            vec[sitetype].push_back(alps::evaluate<double>(matrix_symbolic[i][j]));
        }
    }
  }

  diagonal_matrix_element[name]=vec;
  return true;
}



template <class G, class StateType>
void QMCRun<G,StateType>::initialize_site_states()
{
  maximum_number_of_states=0;
  maximum_sitetype = 0;

  for (site_iterator it=this->sites().first; it!=this->sites().second;++it) {
    unsigned int sitetype=super_type::site_type(*it);
    if (sitetype > maximum_sitetype)
      maximum_sitetype = sitetype;
    site_site_type.push_back(sitetype);
    std::map<int,int>::const_iterator found = number_states_for_site_type_.find(sitetype);
    if (found != number_states_for_site_type_.end())
      site_number_of_states.push_back(found->second);
    else {
      // get site basis
      int num_states = this->model().basis().site_basis(sitetype).num_states();

      if (num_states<=0 || num_states>std::numeric_limits<state_type>::max())
      throw std::invalid_argument("QMC local basis exceeds the state representation");
    number_states_for_site_type_[sitetype]=num_states;
      if (num_states>maximum_number_of_states)
        maximum_number_of_states = num_states;
      site_number_of_states.push_back(num_states);
    }
  }

// #error get rid of the is_charge_model_ and is_spin_model_ below and replace it by the general measurements
  is_charge_model_ = build_diagonal_operator("n");
  is_spin_model_ = build_diagonal_operator("Sz");

  // set up the matrices for the general measurements
  std::map<std::string,std::string>::iterator it = this-> average_expressions.begin();
  while ( it != this-> average_expressions.end() ) {
    std::map<std::string,std::string>::iterator next = it;
    ++next;
    if (!build_diagonal_operator(it->second)) {
      std::clog << "Will not measure \"" << it->first << "\" since it is off-diagonal or not a site operator\n";
      this-> average_expressions.erase(it);
    }
    it = next;
  }

  it = this->local_expressions.begin();
  while ( it != this-> local_expressions.end() ) {
    std::map<std::string,std::string>::iterator next = it;
    ++next;
    if (!build_diagonal_operator(it->second)) {
      std::clog << "Will not measure \"" << it->first << "\" since it is off-diagonal or not a site operator\n";
      this-> local_expressions.erase(it);
    }
    it = next;
  }

  std::map<std::string,std::pair<std::string,std::string> >::iterator itp = this-> correlation_expressions.begin();
  while ( itp != this-> correlation_expressions.end() ) {
    std::map<std::string,std::pair<std::string,std::string> >::iterator next = itp;
    ++next;
    if (!build_diagonal_operator(itp->second.first) || !build_diagonal_operator(itp->second.second)) {
      std::clog << "Will not measure \"" << itp->first << "\" since it is off-diagonal\n";
      this-> correlation_expressions.erase(itp);
    }
    itp = next;
  }

  itp = this-> structurefactor_expressions.begin();
  while ( itp != this-> structurefactor_expressions.end() ) {
    std::map<std::string,std::pair<std::string,std::string> >::iterator next = itp;
    ++next;
    if (!build_diagonal_operator(itp->second.first) || !build_diagonal_operator(itp->second.second)) {
      std::clog << "Will not measure \"" << itp->first << "\" since it is off-diagonal\n";
      this-> structurefactor_expressions.erase(itp);
    }
    itp = next;
  }

}


template <class G, class StateType>
void QMCRun<G,StateType>::create_common_observables()
{
  local.resize(this->num_sites());
  if (is_signed_)
    this->add_measurement("Sign",1,false);

  // Initialize this->measurements
  this->add_measurement("Energy");
  this->add_measurement("Energy Density");

  if (is_charge_model_) {
    measure_local_compressibility_=this->parms.value_or_default(
               "MEASURE[Local Compressibility]",false);
    measure_site_compressibility_=this->parms.value_or_default(
               "MEASURE[Site Compressibility]",false);
    measure_local_density_=this->parms.value_or_default(
                "MEASURE[Local Density]",false)
          || measure_local_compressibility_ || measure_site_compressibility_;
  }

   if (is_spin_model_) {
     measure_local_magnetization_=this->parms.value_or_default("MEASURE[Local Magnetization]",false);
  }


  std::vector<std::string> sitelabels;
  if (measure_local_density_ || measure_local_compressibility_ || measure_local_magnetization_)
    sitelabels =  this->site_labels();

  std::vector<std::string> sitetypelabels;
  if (measure_site_type_density_)
    for (std::size_t i=0; i< num_site_types_;++i)
      sitetypelabels.push_back(boost::lexical_cast<std::string>(i));

  std::vector<std::string> bondtypelabels;
  std::vector<std::string> bondtypecorrlabels;
  if (measure_bond_type_stiffness_) {
    for (std::size_t i=0; i< num_bond_types_;++i) {
      bondtypelabels.push_back(boost::lexical_cast<std::string>(i));
      for (std::size_t j=0; j< num_bond_types_;++j)
        bondtypecorrlabels.push_back(boost::lexical_cast<std::string>(i) + " -- "
                                   + boost::lexical_cast<std::string>(j));
    }
  }

  std::vector<std::string> corrlabels;
  if (measure_correlations_ || measure_green_function_ || !this->correlation_expressions.empty()) {
    if (measurement_origin_) {
      corrlabels.clear();
      for (unsigned int j=0;j<this->num_sites();++j)
        corrlabels.push_back(this->coordinate_string(measurement_origin_.get())+" -- " + this->coordinate_string(j));
    }
    else {
      corrlabels = this->distance_labels();
      distance_mult = this->distance_multiplicities();
    }
  }


  std::vector<std::string> momentalabels;
  if(measure_structure_factor_ || !this->structurefactor_expressions.empty())
    momentalabels = this->momenta_labels();

  if (winding_dimension()) this->add_measurement("Stiffness");
  if (measure_bond_type_stiffness_) {
    this->add_measurement("Bond Type Stiffness",bondtypelabels);
    this->add_measurement("Bond Type Stiffness Correlations",bondtypecorrlabels);
  }
  if (is_spin_model_) {
    this->add_measurement("Magnetization");
    this->add_measurement("Magnetization Density");
    this->add_measurement("|Magnetization|");
    this->add_measurement("|Magnetization Density|");
    this->add_measurement("Magnetization^2");
    this->add_measurement("Magnetization Density^2");
    this->add_measurement("Magnetization^4");
    this->add_measurement("Magnetization Density^4");
    this->add_measurement("Susceptibility");
    if (this->is_bipartite()) {
      this->add_measurement("Staggered Magnetization");
      this->add_measurement("Staggered Magnetization Density");
      this->add_measurement("|Staggered Magnetization|");
      this->add_measurement("|Staggered Magnetization Density|");
      this->add_measurement("Staggered Magnetization^2");
      this->add_measurement("Staggered Magnetization Density^2");
      this->add_measurement("Staggered Magnetization^4");
      this->add_measurement("Staggered Magnetization Density^4");
    }
    if(measure_correlations_)
      this->add_measurement("Spin Correlations",corrlabels);
    if(measure_structure_factor_)
      this->add_measurement("Spin Structure Factor",momentalabels);
    if (measure_local_magnetization_)
      this->add_measurement("Local Magnetization",sitelabels);
  }

if(measure_CBS_order_)
{
  this->add_measurement("CBS Order");
  this->add_measurement("CBS Order^2");
  this->add_measurement("CBS Order^4");
}

  if (is_charge_model_) {
    this->add_measurement("Centered Density Moments",3,false);
    this->add_measurement("Density");
    this->add_measurement("Density^2");
    if (this->is_bipartite()) {
      this->add_measurement("Checkerboard Charge Order");
      this->add_measurement("Checkerboard Charge Order^2");
      this->add_measurement("Checkerboard Charge Order^4");
    }
    if (measure_local_density_)
      this->add_measurement("Local Density",sitelabels);
    if (measure_site_type_density_)
      this->add_measurement("Site Type Density",sitetypelabels);
    if (measure_site_compressibility_)
      this->add_measurement("Integrated Local Density Correlations",sitelabels);
    if (measure_local_compressibility_)
      this->add_measurement("Local Density * Global Density",sitelabels);
    if(measure_correlations_)
      this->add_measurement("Density Correlations",corrlabels);
    if(measure_structure_factor_)
      this->add_measurement("Density Structure Factor",momentalabels);
  }

  if (measure_green_function_)
    this->add_measurement("Green's Function",corrlabels);

  // set up the matrices for the general measurements
  std::map<std::string,std::string>::iterator it = this-> average_expressions.begin();
  while ( it != this-> average_expressions.end() ) {
    std::map<std::string,std::string>::iterator next = it;
    ++next;
    if (this->measurements.count(it->first)) {
      std::clog << "Will not measure custom average measurement \"" << it->first << "\" since a hard-coded measurement already exists\n";
      this-> average_expressions.erase(it);
    }
    else
      this->add_measurement(it->first);
    it = next;
  }

  it = this->local_expressions.begin();
  while ( it != this-> local_expressions.end() ) {
    std::map<std::string,std::string>::iterator next = it;
    ++next;
    if (this->measurements.count(it->first)) {
      std::clog << "Will not measure custom local measurement \"" << it->first << "\" since a hard-coded measurement already exists\n";
      this-> local_expressions.erase(it);
    }
    else {
      if (sitelabels.empty())
        sitelabels = this->site_labels();
      this->add_measurement(it->first,sitelabels);
    }
    it = next;
  }

  std::map<std::string,std::pair<std::string,std::string> >::iterator itp = this-> correlation_expressions.begin();
  while ( itp != this-> correlation_expressions.end() ) {
   std::map<std::string,std::pair<std::string,std::string> >::iterator next = itp;
    ++next;
    if (this->measurements.count(itp->first)) {
      std::clog << "Will not measure custom correlation measurements \"" << itp->first << "\" since a hard-coded measurement already exists\n";
      this-> correlation_expressions.erase(itp);
    }
    else
      this->add_measurement(itp->first,corrlabels);
    itp = next;
  }

  itp = this-> structurefactor_expressions.begin();
  while ( itp != this-> structurefactor_expressions.end() ) {
   std::map<std::string,std::pair<std::string,std::string> >::iterator next = itp;
    ++next;
    if (this->measurements.count(itp->first)) {
      std::clog << "Will not measure custom correlation measurements \"" << itp->first << "\" since a hard-coded measurement already exists\n";
      this-> structurefactor_expressions.erase(itp);
    }
    else
      this->add_measurement(itp->first,momentalabels);
    itp = next;
  }
}


template <class G, class StateType>
bool QMCRun<G,StateType>::do_common_measurements(double sign, const std::vector<state_type>& state, const std::valarray<double>& localint)
{
  std::vector<std::vector<double> > const& matrix_element_Sz = diagonal_matrix_element["Sz"];
  std::vector<std::vector<double> > const& matrix_element_n = diagonal_matrix_element["n"];

  double NbSites=this->num_sites();

  auto total=[&](auto const& table) {
    double value=0;
    for (size_t i=0;i<this->num_sites();++i) value+=table[this->site_type(i)][state[i]];
    return value;
  };
  if ((restricted_magnetization && (!is_spin_model_ || *restricted_magnetization!=total(matrix_element_Sz))) ||
      (restricted_particle_number && (!is_charge_model_ || *restricted_particle_number!=total(matrix_element_n))))
    return false;
  if (is_signed_)
    this->record("Sign",sign,sign);

  if (is_spin_model_) {
    double sz=0.;
    double ssz=0.;
    for (unsigned int i=0;i<this->num_sites();++i) {
      double x = local[i] = matrix_element_Sz[this->site_type(i)][state[i]];
      sz+=x;
      if (this->is_bipartite())
        ssz += this->parity(i)*x;
    }


    this->record("Magnetization",sz*sign,sign);
    this->record("Magnetization Density",sz/NbSites*sign,sign);
    this->record("|Magnetization|",std::abs(sz)*sign,sign);
    this->record("|Magnetization Density|",std::abs(sz)/NbSites*sign,sign);
    this->record("Magnetization^2",sz*sz*sign,sign);
    this->record("Magnetization Density^2",sz*sz/NbSites/NbSites*sign,sign);
    this->record("Magnetization^4",sz*sz*sz*sz*sign,sign);
    this->record("Magnetization Density^4",sz*sz*sz*sz/NbSites/NbSites/NbSites/NbSites*sign,sign);
    this->record("Susceptibility",sz*sz*beta/NbSites*sign,sign);
    if (this->is_bipartite()) {
      this->record("Staggered Magnetization",ssz*sign,sign);
      this->record("Staggered Magnetization Density",ssz/NbSites*sign,sign);
      this->record("|Staggered Magnetization|",std::abs(ssz)*sign,sign);
      this->record("|Staggered Magnetization Density|",std::abs(ssz)/NbSites*sign,sign);
      this->record("Staggered Magnetization^2",ssz*ssz*sign,sign);
      this->record("Staggered Magnetization Density^2",ssz*ssz/NbSites/NbSites*sign,sign);
      this->record("Staggered Magnetization^4",ssz*ssz*ssz*ssz*sign,sign);
      this->record("Staggered Magnetization Density^4",ssz*ssz*ssz*ssz/NbSites/NbSites/NbSites/NbSites*sign,sign);
    }
    if (measure_local_magnetization_) {
      std::valarray<double> local_sz(local);
      local_sz *= sign;
      this->record("Local Magnetization",local_sz,sign);
    }
  }

  if (is_charge_model_) {
    if(measure_site_type_density_) {
      site_type_density_.resize(num_site_types_);
      for (int i=0;i<num_site_types_;++i)
        site_type_density_[i]=0.;
    }
    double n=0.;
    double sn=0.;
    for (unsigned int i=0;i<this->num_sites();++i) {
      double nloc = local[i] = matrix_element_n[this->site_type(i)][state[i]];
      n+=nloc;
      if(measure_site_type_density_)
        site_type_density_[this->site_type(i)] += nloc;
      if (this->is_bipartite())
        sn += this->parity(i)*nloc;

    }


    if (this->is_bipartite()) {
      this->record("Checkerboard Charge Order",sn/NbSites*sign,sign);
      this->record("Checkerboard Charge Order^2",sn*sn/NbSites/NbSites*sign,sign);
      this->record("Checkerboard Charge Order^4",sn*sn*sn*sn/NbSites/NbSites/NbSites/NbSites*sign,sign);
    }
    if(measure_site_type_density_) {
      site_type_density_ *= sign/this->num_sites();
      this->record("Site Type Density",site_type_density_,sign);
    }
    if (measure_local_density_) {
      std::valarray<double> local_density(local);
      local_density *= sign;
      this->record("Local Density",local_density,sign);
      if (measure_local_compressibility_) {
        local_density *= n/NbSites;
        this->record("Local Density * Global Density",local_density,sign);
      }
      if (measure_site_compressibility_) {
        BOOST_ASSERT (local_density.size() == localint.size());
        local_density=localint*sign;
        this->record("Integrated Local Density Correlations",local_density,sign);
      }
    }
    this->record("Density",n/NbSites*sign,sign);
    this->record("Density^2",n*n/NbSites*sign,sign);
  }

  // custom average measurements

  typedef std::pair<std::string,std::string> string_pair_type;
  BOOST_FOREACH(string_pair_type const& x, this->average_expressions) {
    std::vector<std::vector<double> > const& matrix_element = diagonal_matrix_element[x.second];
    double ave=0.;
    for (unsigned int i=0;i<this->num_sites();++i)
      ave += matrix_element[this->site_type(i)][state[i]];
    this->record(x.first,ave/NbSites*sign,sign);
  }

  // custom local measurements
  BOOST_FOREACH(string_pair_type const& x, this->local_expressions) {
    std::vector<std::vector<double> > const& matrix_element = diagonal_matrix_element[x.second];
    std::valarray<double> local(this->num_sites());
    for (unsigned int i=0;i<this->num_sites();++i)
      local[i] = matrix_element[this->site_type(i)][state[i]];
    local *= sign;
    this->record(x.first,local,sign);
  }

  // custom correlation measurements
  typedef std::pair<std::string,std::pair<std::string,std::string> > string_pair_pair_type;
  BOOST_FOREACH(string_pair_pair_type const& x, this->correlation_expressions) {
    std::vector<std::vector<double> > const& matrix_element_a = diagonal_matrix_element[x.second.first];
    std::vector<std::vector<double> > const& matrix_element_b = diagonal_matrix_element[x.second.second];
    std::valarray<double> corr(measurement_origin_ ? this->num_sites() : this->num_distances());
    std::vector<double> local_a(this->num_sites());
    std::vector<double> local_b(this->num_sites());
    for (int i=0;i<this->num_sites();++i) {
      local_a[i] = matrix_element_a[this->site_type(i)][state[i]];
      local_b[i] = matrix_element_b[this->site_type(i)][state[i]];
    }
    corr=0.;
    if (measurement_origin_) {
      for (size_t i=0;i<this->num_sites();++i) corr[i]=local_a[*measurement_origin_]*local_b[i]*sign;
    } else {
      for (size_t i=0;i<this->num_sites();++i)
        for (size_t j=0;j<this->num_sites();++j)
          corr[this->distance(i,j)]+=local_a[i]*local_b[j];
      for (size_t i=0;i<corr.size();++i) corr[i]*=sign/distance_mult[i];
    }
    this->record(x.first,corr,sign);
  }

  // Correlation measurements

  if(measure_correlations_) {
    std::valarray<double> corr;
    if(measurement_origin_) {
      corr.resize(this->num_sites());
      int i1=measurement_origin_.get();
      for (int i2=0;i2<this->num_sites();++i2)
        corr[i2] += local[i1]*local[i2];
    }
    else {
      corr.resize(this->num_distances());
      corr=0.;
      for (int i1=0;i1<this->num_sites();++i1)
        for (int i2=0;i2<this->num_sites();++i2) {
          typename super_type::size_type d = this->distance(i1,i2);
          BOOST_ASSERT(d>=0 && d< corr.size());
          corr[this->distance(i1,i2)] += local[i1]*local[i2];
        }
      for (int i=0;i<corr.size();++i)
        corr[i]/=distance_mult[i];
    }
    corr *= sign;
    this->record(is_charge_model_ ? "Density Correlations" : "Spin Correlations",corr,sign);
  }

  // Structure factor measurements

  // custom structure factor measurements
  BOOST_FOREACH(string_pair_pair_type const& x, this->structurefactor_expressions) {
    std::vector<std::vector<double> > const& matrix_element_a = diagonal_matrix_element[x.second.first];
    std::vector<std::vector<double> > const& matrix_element_b = diagonal_matrix_element[x.second.second];
    std::vector<double> str;
    str.clear();
    std::vector<double> local_a(this->num_sites());
    std::vector<double> local_b(this->num_sites());
    for (int i=0;i<this->num_sites();++i) {
      local_a[i] = matrix_element_a[this->site_type(i)][state[i]];
      local_b[i] = matrix_element_b[this->site_type(i)][state[i]];
    }
    for (typename super_type:: momentum_iterator mit=this->momenta().first; mit != this->momenta().second; ++mit) {
      std::complex<double> vala, valb;
      for (typename super_type::site_iterator sit=this->sites().first; sit!=this->sites().second;++sit) {
        double phase = alps::numeric::scalar_product(this->momentum(*mit), super_type::coordinate(*sit));
        std::complex<double> cphase(std::cos(phase), std::sin(phase));
        vala += local_a[*sit] * cphase;
        valb += local_b[*sit] * cphase;
      }
      str.push_back(std::real(std::conj(vala)*valb));
    }
    std::valarray<double> strv(str.size());
    for (int i=0;i<str.size();++i)
      strv[i]=str[i]*sign/this->num_sites();
    this->record(x.first,strv,sign);
  }

  if(measure_structure_factor_) {
    std::vector<double> str;
    for (typename super_type:: momentum_iterator mit=this->momenta().first; mit != this->momenta().second; ++mit) {
      std::complex<double> val;
      for (typename super_type::site_iterator sit=this->sites().first; sit!=this->sites().second;++sit) {
        double phase = alps::numeric::scalar_product(this->momentum(*mit), super_type::coordinate(*sit));
        val += local[*sit] * std::complex<double>(std::cos(phase), std::sin(phase));
      }
      str.push_back(std::norm(val));
    }
    std::valarray<double> strv(str.size());
    for (int i=0;i<str.size();++i)
      strv[i]=str[i]*sign/this->num_sites();
    this->record(is_charge_model_ ? "Density Structure Factor" : "Spin Structure Factor",strv,sign);
  }
  if (measure_CBS_order_)
  {
        std::complex<double> val;
        for (typename super_type::site_iterator sit=this->sites().first; sit!=this->sites().second;++sit)
        {
            double phase = boost::math::constants::pi<double>() * (super_type::coordinate(*sit)[0] + super_type::coordinate(*sit)[1]);
            val += local[*sit] * (std::complex<double> ( std::cos(phase), std::sin(phase) ));
        }
        double CBS_order = (std::abs(val) * std::abs(val)) / this->num_sites();
        this->record("CBS Order",CBS_order*sign,sign);
        this->record("CBS Order^2",CBS_order * CBS_order * sign,sign);
        this->record("CBS Order^4",CBS_order * CBS_order * CBS_order * CBS_order * sign,sign);
  }

  return true;
}

#endif
