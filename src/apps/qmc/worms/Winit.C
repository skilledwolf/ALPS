/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2001-2005 by Matthias Troyer <troyer@comp-phys.org>,
*                            Simon Trebst <trebst@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include "WRun.h"
#include <iomanip>
//#include <alps/osiris/std/set.h>

void WRun::print_copyright(std::ostream& out)
{
  out << "Worm algorithm quantum Monte Carlo simulation v2.0\n"
      << "  copyright (c) 1997-2007 by Simon Trebst <trebst@comp-phys.org>\n"
      << "                          and Matthias Troyer <troyer@comp-phys.org>\n"
      << " for details see the publication:\n"     
      << " A.F. Albuquerque et al., J. of Magn. and Magn. Materials 310, 1187 (2007).\n\n";
}



//- Constructor/setup run -------------------------------------------------

WRun::WRun(alps::params const& p,size_t bins,size_t chain)
  : QMCRun<>(p,bins,chain,p.exists("NUMBER_OF_PARTICLES")),
    chain_(chain),
    canonical(parms.defined("NUMBER_OF_PARTICLES")),
    adjust_parameter(parms.defined("ADJUST") ? static_cast<std::string>(parms["ADJUST"]) : std::string("mu")),
    worm_head(2,wormhead_type(graph(),kinks)),
    stat(4),
    eta(beta*double(parms.value_or_default("eta", 1.))),
    thermal_sweeps(parms.required_value("THERMALIZATION")),
    skip_measurements(parms.value_or_default("SKIP",1)),
    have_worm(false),
    chain_kappa(parms.value_or_default("CHAIN_KAPPA",false)),
    worms_per_kink(parms.value_or_default("WORMS_PER_KINK",1)),
    log_numeric_limits_double(std::log(std::numeric_limits<double>::max())),
    nonlocal(parms.value_or_default("NONLOCAL",true)),
    use_1D_stiffness(parms.value_or_default("USE_1D_STIFFNESS",false)),  //@#$br
    chain_number(num_sites()),
    bond_type(alps::get_or_default(alps::bond_type_t(),graph(),0))
{

  if (measure_green_function_ || parms.value_or_default("MEASURE[Site Compressibility]",false) || measure_bond_type_stiffness_)
    throw std::invalid_argument("Worm Green function, site compressibility and bond-type stiffness estimators are not implemented");
  if (use_1D_stiffness && winding_dimension()!=1)
    throw std::invalid_argument("USE_1D_STIFFNESS requires a one-dimensional lattice");
  if(!std::isfinite(eta) || eta<=0)
    boost::throw_exception(std::out_of_range("negative eta is illegal"));
  if(thermal_sweeps<0)
    boost::throw_exception( std::out_of_range("negative thermalization is illegal"));
  
  // kinks
  {
    int vol = num_sites();
    kinks.resize(vol);
    initial_state_.resize(vol,state_type());
  }
  
  // model
  initialize_site_states();
  if (canonical && (!is_charge_model_ || !parms.defined("CORRECTION")))
    throw std::invalid_argument("Particle-number adjustment requires a charge model and CORRECTION");
  original_bond_type.resize(num_bonds());
  for (auto [it,end]=bonds();it!=end;++it) original_bond_type[index(*it)]=bond_type[*it];
  initialize_hamiltonian();
  if (canonical) {
    if (!std::isfinite(double(parms["CORRECTION"])) || double(parms["CORRECTION"])<=0)
      throw std::invalid_argument("CORRECTION must be finite and positive");
    double lower=0., upper=0.;
    for (auto [it,end]=sites();it!=end;++it) {
      auto const& n=QMCRun<>::diagonal_matrix_element.at("n")[site_type(*it)];
      lower+=*std::min_element(n.begin(),n.end());
      upper+=*std::max_element(n.begin(),n.end());
    }
    double target=double(parms["NUMBER_OF_PARTICLES"]);
    if (target<lower || target>upper)
      throw std::invalid_argument("NUMBER_OF_PARTICLES is outside the model's allowed range");
    requested_site_matrix=site_matrix;
  }
  // print_hamiltonian();

  // subintervals for non-local interactions
  subinterval.clear();
  subinterval_valid = false;
  current_head_num = 0;

  num_kinks=num_sites();
  create_observables();
  start();
}


namespace {
template<class State,class Transfer> void checkpoint_fields(State& s,Transfer field) {
  field("steps",s.steps);
  field("measurement_cursor",s.measurements_done);
  field("num_kinks",s.num_kinks);
  field("worms_per_update",s.worms_per_update);
  field("sign",s.Sign);
  field("last_id",s.last_id_);
  field("initial_state",s.initial_state_);
  field("corrections_upwards",s.corrections_upwards);
  field("corrections_downwards",s.corrections_downwards);
  field("preadjustment_done",s.preadjustment_done);
  field("adjustment_done",s.adjustment_done);
  field("adjustment_samples",s.nob);
}
}

std::vector<double> WRun::hamiltonian_state(bool include_sites) const {
  std::vector<double> values=hopping_matrix;
  if (include_sites) for (auto const& row:site_matrix) values.insert(values.end(),row.begin(),row.end());
  for (auto const& [type,matrix]:diagonal_matrix_element)
    values.insert(values.end(),matrix.data(),matrix.data()+matrix.num_elements());
  for (auto const& row:matrix_element_raise_) values.insert(values.end(),row.begin(),row.end());
  for (auto const& row:matrix_element_lower_) values.insert(values.end(),row.begin(),row.end());
  for (auto [it,end]=bonds();it!=end;++it) {
    values.push_back(source(*it)); values.push_back(target(*it)); values.push_back(bond_type[*it]);
  }
  return values;
}

void WRun::save(alps::hdf5::archive& ar) const {
  QMCRun<>::save(ar);
  ar["checkpoint/version"] << uint64_t(1);
  ar["checkpoint/chain"] << uint64_t(chain_);
  ar["checkpoint/hamiltonian"] << hamiltonian_state();
  if (canonical) ar["checkpoint/adjusted_parameter"] << double(parms[adjust_parameter]);
  checkpoint_fields(static_cast<worm_state const&>(*this),[&](char const* key,auto const& value){
    ar[std::string("checkpoint/")+key] << value;
  });
  std::vector<std::array<uint32_t,3>> index;
  std::vector<double> times;
  for (size_t site=0;site<kinks.size();++site) for (auto const& kink:kinks[site]) {
    index.push_back({uint32_t(site),kink.state(),kink.id()}); times.push_back(kink.time());
  }
  ar["checkpoint/kinks/index"] << index;
  ar["checkpoint/kinks/time"] << times;
}

void WRun::load(alps::hdf5::archive& ar) {
  alps::params saved;
  uint64_t version,chain;
  ar["/parameters"] >> saved;
  ar["checkpoint/version"] >> version;
  ar["checkpoint/chain"] >> chain;
  if (version!=1 || chain!=chain_ || checkpoint_parameters(saved)!=checkpoint_parameters(parameters))
    throw std::invalid_argument("Worm checkpoint model or version mismatch");
  WRun staged(parameters,bins_,chain_);
  if (canonical) {
    double value;
    ar["checkpoint/adjusted_parameter"] >> value;
    if (!std::isfinite(value)) throw std::invalid_argument("Invalid adjusted worm parameter");
    staged.set_adjusted_parameter(value);
  }
  std::vector<double> hamiltonian;
  ar["checkpoint/hamiltonian"] >> hamiltonian;
  if (hamiltonian!=staged.hamiltonian_state())
    throw std::invalid_argument("Worm checkpoint lattice or Hamiltonian changed");
  checkpoint_fields(static_cast<worm_state&>(staged),[&](char const* key,auto& value){
    ar[std::string("checkpoint/")+key] >> value;
  });
  auto production=staged.steps>thermal_sweeps ? staged.steps-thermal_sweeps : 0;
  if (production>parameters["SWEEPS"].as<uint64_t>() ||
      staged.measurements_done<=0 || staged.measurements_done>skip_measurements ||
      !std::isfinite(staged.worms_per_update) || staged.worms_per_update<0 || staged.worms_per_update>INT_MAX ||
      (staged.Sign!=1. && staged.Sign!=-1.) || (!is_signed_ && staged.Sign!=1.) ||
      staged.corrections_upwards<0 || staged.corrections_downwards<0 ||
      staged.nob.size()>=60 ||
      (staged.adjustment_done && !staged.preadjustment_done) ||
      ((!staged.preadjustment_done || staged.adjustment_done) && !staged.nob.empty()) ||
      (canonical && staged.steps>=thermal_sweeps && !staged.adjustment_done) ||
      staged.initial_state_.size()!=num_sites())
    throw std::invalid_argument("Invalid worm checkpoint progress");
  std::vector<std::array<uint32_t,3>> index;
  std::vector<double> times;
  ar["checkpoint/kinks/index"] >> index;
  ar["checkpoint/kinks/time"] >> times;
  if (index.size()!=times.size() || index.size()%2 || staged.num_kinks!=num_sites()+index.size()/2)
    throw std::invalid_argument("Invalid worm checkpoint kink count");
  for (size_t i=0;i<index.size();++i) {
    auto [site,state,id]=index[i];
    if (site>=num_sites() || state>=site_number_of_states[site] || id>staged.last_id_ ||
        !std::isfinite(times[i]) || times[i]<0 || times[i]>=1 ||
        (!staged.kinks[site].empty() && double(staged.kinks[site].back().time())>=times[i]))
      throw std::invalid_argument("Invalid worm checkpoint kink");
    kink_type kink(times[i],state); kink.set_id(id); staged.kinks[site].push_back(kink);
  }
  struct jump { size_t site; double time; int change; };
  std::map<unsigned,std::vector<jump>> jumps;
  for (size_t site=0;site<num_sites();++site) {
    if (staged.initial_state_[site]>=site_number_of_states[site])
      throw std::invalid_argument("Invalid worm initial state");
    if (staged.kinks[site].empty()) continue;
    int before=staged.kinks[site].back().state();
    for (auto const& kink:staged.kinks[site]) {
      int change=int(kink.state())-before;
      if (std::abs(change)!=1) throw std::invalid_argument("Invalid worm worldline jump");
      jumps[kink.id()].push_back({site,double(kink.time()),change}); before=kink.state();
    }
  }
  for (auto const& [id,pair]:jumps) {
    if (pair.size()!=2 || pair[0].time!=pair[1].time || pair[0].change!=-pair[1].change || pair[0].site==pair[1].site)
      throw std::invalid_argument("Unpaired worm worldline jump");
    bool adjacent=false;
    for (unsigned n=0;n<num_neighbors(pair[0].site);++n)
      adjacent=adjacent || neighbor(pair[0].site,n)==pair[1].site;
    if (!adjacent) throw std::invalid_argument("Worm jump does not follow a lattice bond");
  }
  validate_measurements(ar,production/skip_measurements);
  staged.QMCRun<>::load(ar);
  // Publish only after all physical state, native accumulators and RNG validate.
  static_cast<worm_state&>(*this)=std::move(static_cast<worm_state&>(staged));
  measurements=std::move(staged.measurements); random=std::move(staged.random);
  parms=std::move(staged.parms);
  adjustment_energy_shift=staged.adjustment_energy_shift;
  density_reference_=staged.density_reference_;
  site_matrix=std::move(staged.site_matrix); hopping_matrix=std::move(staged.hopping_matrix);
  diagonal_matrix_element=std::move(staged.diagonal_matrix_element);
  matrix_element_raise_=std::move(staged.matrix_element_raise_);
  matrix_element_lower_=std::move(staged.matrix_element_lower_);
  subinterval_valid=false; have_worm=false;
}
