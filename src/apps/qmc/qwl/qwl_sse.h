/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2004-2005 by Stefan Wessel <wessel@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#ifndef ALPS_QWL_SSE_H
#define ALPS_QWL_SSE_H

#include "qwl_histogram.h"
#include <alps/mcbase.hpp>
#include <alps/model.h>
#include <alps/ngs/make_deprecated_parameters.hpp>
#include <alps/alea/checkpoint.hpp>
#include <alps/alea/hdf5.hpp>
#include <alps/hdf5/stdarray.hpp>
#include <vector>
#include <algorithm>
#include <string>
#include <cmath>

using namespace std;
using namespace alps;

template<typename T> T sqr(T x) {return x*x;}


struct qwl_state {
  uint64_t sweeps=0, production_sweeps=0, thistime=0, block_sweeps=0, block_sweeps_total=0;
  unsigned int non0=0, norder=0, logf_step=1;
  double logf=0, minimum_histogram=0, flatness_treshold=0;
  bool doing_multicanonical=false, upwalker=true, all_done=false;
  histogram<double> g, histo, histoup, nmeasurements, uniform_structure_factor,
                    staggered_structure_factor, transition_prob;
  struct vertextype {
    uint32_t op=0, bond=0, vs=0, vn[4]{}, vvisited[4]{};
    int vv[4]{};
  };
  std::vector<vertextype> operator_string;
  std::vector<uint32_t> state;
};

class QWL_SSE_Simulation : public alps::mcbase, private alps::graph_helper<>,
                           private alps::model_helper<>, private qwl_state {
public:
  QWL_SSE_Simulation(alps::params const&, size_t bins=128, size_t chain=0);
  void save(alps::hdf5::archive&) const override;
  void load(alps::hdf5::archive&) override;
  void update() override;
  void measure() override {} // Sampling occurs inside the extended-ensemble update.
  double fraction_completed() const override { return all_done ? 1. : 0.; }
  uint64_t completed_sweeps() const { return sweeps; }
  size_t number_of_sites() const { return num_sites(); }
  bool bipartite() const { return is_bipartite(); }
  static alps::params checkpoint_parameters(alps::params p) { p.erase("SWEEPS"); return p; }
private:
  using statetype=uint32_t;
  using bond_state_type=uint32_t;
  static alps::Parameters graph_parameters(alps::params const& p) {
    alps::Disorder::seed(p.value_or<uint32_t>("DISORDER_SEED",uint32_t(p.value_or("SEED",42))));
    return alps::make_deprecated_parameters(p);
  }
  alps::Parameters parms;
  size_t bins_, chain_;
  unsigned int L, norder_min, norder_max, logf_steps_total;
  double beta=1, offset=0, initial_logf=0;
  bool includeLfactor, use_zhou_bhatt, measure_magnetics;
  std::vector<int> statev;
  std::vector<uint32_t> staten;
  double random_01() const { return random(); }
  int random_int(int a,int b) { return a+int((b-a+1)*random()); }
  void record(std::string const& name,std::valarray<double> const& value) {
    record(name,std::vector<double>(std::begin(value),std::end(value)));
  }
  template<class T> void record(std::string const& name,T const& value) {
    auto sample=alps::alea::make_adapter(value);
    if (name.compare(0,5,"Time ")==0) {
      *measurement(name)<<sample;
      *measurement<alps::alea::autocorr_acc<double>>("Autocorrelation: "+name)<<sample;
    }
    else *measurement<alps::alea::mean_acc<double>>(name)<<sample;
  }
  std::vector<std::array<uint64_t,3>> topology() const {
    std::vector<std::array<uint64_t,3>> value;
    for (auto [it,end]=bonds();it!=end;++it)
      value.push_back({uint64_t(source(*it)),uint64_t(target(*it)),uint64_t(inhomogeneous_bond_type(*it))});
    return value;
  }
  vector <boost::uint32_t> is_ferromagnetic;
  vector <double> matrix_factor;

  void init_measurements();
  void init_tables();
  void diagonal_update();
  void diagonal_update_multicanonical();
  void vertexbuild();
  void loopupdate();
  void measure_magnetism();
  void store_histograms(string);
  boost::uint32_t exit_leg(vertextype& , boost::uint32_t);
  double matrix_element(const int, const int, const int);

  statetype get_state(bond_state_type vs, int n) const {
    return ((vs>>n) & 1);
  }

  void set_statevector_in(bond_state_type& vs, statetype s0, statetype s1) {
    vs=s0 | (s1<<1);
  }

  void set_statevector_out(bond_state_type& vs, statetype s2, statetype s3) {
    vs=(vs | (s2<<2) ) | (s3<<3);
  }

  double decode_spin(statetype s) const {
    return ( (s) ? 0.5 : -0.5);
  }

  statetype random_state() const {
    return (random_01()>0.5 ? 1 : 0);
  }

  void transform(bond_state_type& vs, int n)  {
    vs^=(1<<n);
  }

};

QWL_SSE_Simulation::QWL_SSE_Simulation(alps::params const& p,size_t bins,size_t chain)
  : mcbase(p,chain), graph_helper<>(graph_parameters(p)),
    model_helper<>(static_cast<graph_helper<> const&>(*this),alps::make_deprecated_parameters(p)),
    parms(alps::make_deprecated_parameters(p)), bins_(bins), chain_(chain)
{
  if (!num_sites() || !num_bonds() || num_sites()>INT_MAX || num_bonds()>INT_MAX || bins<2 || bins%2)
    throw std::invalid_argument("QWL requires a nonempty lattice and even batch capacity >= 2");
  state.resize(num_sites());
  statev.resize(num_sites());
  staten.resize(num_sites());
  for (site_iterator it=sites().first;it!=sites().second;++it)
    state[*it]=random_state();
  logf_steps_total=parms.value_or_default("NUMBER_OF_WANG_LANDAU_STEPS",16);
  norder_min=parms.value_or_default("EXPANSION_ORDER_MINIMUM",0);
  norder_max=parms.value_or_default("EXPANSION_ORDER_MAXIMUM",parms.value_or_default("CUTOFF",500));
  if (norder_min>norder_max || norder_max>=INT_MAX || !logf_steps_total || logf_steps_total>64)
    throw std::invalid_argument("Invalid QWL expansion window or Wang-Landau step count");
  int g_left=norder_min;
  int g_size=norder_max-g_left+1;
  if (g_left<0)
    boost::throw_exception(std::runtime_error("EXPANSION_ORDER_MINIMUM must not be negative"));
  g.resize(g_size,g_left);
  histo.resize(g_size,g_left);
  histoup.resize(g_size,g_left);
  nmeasurements.resize(g_size,g_left);
  transition_prob.resize(g_size,g_left);
  L=parms.value_or_default("CUTOFF",norder_max);
  if (!L || L>=INT_MAX || L<norder_max)
    boost::throw_exception(std::runtime_error("EXPANSION_ORDER_MAXIMUM must not exceed CUTOFF"));
  operator_string.resize(L);
  measure_magnetics=parms.value_or_default("MEASURE_MAGNETIC_PROPERTIES",1);
  if (measure_magnetics) {
    uniform_structure_factor.resize(g_size,g_left);
    if (is_bipartite()) staggered_structure_factor.resize(g_size,g_left);
  }
  includeLfactor=parms.value_or_default("INCLUDE_COMBINATORICS_FACTORS",1);
  use_zhou_bhatt=parms.value_or_default("USE_ZHOU_BHATT_METHOD",1.);
  if (use_zhou_bhatt) {
    block_sweeps_total=g.size();
    logf=log(double(parms.value_or_default("INITIAL_MODIFICATION_FACTOR",exp(1.))));
    minimum_histogram=1./logf;
    flatness_treshold=parms.value_or_default("FLATNESS_TRESHOLD",1E10);
  }
  else {
    block_sweeps_total=parms.value_or_default("BLOCK_SWEEPS",10000);
    logf=log(double(parms.value_or_default("INITIAL_INCREASE_FACTOR",exp(L*log((double)num_sites())/block_sweeps_total))));
    minimum_histogram=0.;
    flatness_treshold=parms.value_or_default("FLATNESS_TRESHOLD",0.2);
  }
  initial_logf=logf;
  init_tables();
  if (!block_sweeps_total || !std::isfinite(logf) || logf<=0 ||
      !std::isfinite(minimum_histogram) || !std::isfinite(flatness_treshold) || flatness_treshold<=0)
    throw std::invalid_argument("Invalid QWL refinement controls");
  // A window may start above zero; seed a closed diagonal string in its range.
  if (norder_min) {
    auto it=bonds().first;
    while (it!=bonds().second && matrix_factor[inhomogeneous_bond_type(*it)]==0) ++it;
    if (it==bonds().second) throw std::invalid_argument("QWL window requires a nonzero coupling");
    state[source(*it)]=0;
    state[target(*it)]=is_ferromagnetic[inhomogeneous_bond_type(*it)] ? 0 : 1;
    for (unsigned i=0;i<norder_min;++i) { operator_string[i].op=1; operator_string[i].bond=static_cast<uint32_t>(std::distance(bonds().first,it)); }
    non0=norder=norder_min;
  }
  init_measurements();
}

void QWL_SSE_Simulation::save(alps::hdf5::archive& ar) const {
  mcbase::save(ar);
  ar["checkpoint/version"] << uint64_t(1);
  ar["checkpoint/chain"] << uint64_t(chain_);
  ar["checkpoint/bonds"] << topology();
  ar["checkpoint/factors"] << matrix_factor;
  ar["checkpoint/ferromagnetic"] << is_ferromagnetic;
  ar["checkpoint/sweeps"] << sweeps;
  ar["checkpoint/production_sweeps"] << production_sweeps;
  ar["checkpoint/thistime"] << thistime;
  ar["checkpoint/non0"] << non0;
  ar["checkpoint/norder"] << norder;
  ar["checkpoint/logf_step"] << logf_step;
  ar["checkpoint/block_sweeps"] << block_sweeps;
  ar["checkpoint/block_sweeps_total"] << block_sweeps_total;
  ar["checkpoint/logf"] << logf;
  ar["checkpoint/minimum_histogram"] << minimum_histogram;
  ar["checkpoint/flatness_treshold"] << flatness_treshold;
  ar["checkpoint/doing_multicanonical"] << doing_multicanonical;
  ar["checkpoint/upwalker"] << upwalker;
  ar["checkpoint/all_done"] << all_done;
  ar["checkpoint/state"] << state;
  ar["checkpoint/g"] << g;
  ar["checkpoint/histo"] << histo;
  ar["checkpoint/histoup"] << histoup;
  ar["checkpoint/nmeasurements"] << nmeasurements;
  ar["checkpoint/uniform_structure_factor"] << uniform_structure_factor;
  ar["checkpoint/staggered_structure_factor"] << staggered_structure_factor;
  ar["checkpoint/transition_prob"] << transition_prob;
  std::vector<std::array<uint32_t,2>> ops;
  for (auto const& v:operator_string) ops.push_back({v.op,v.bond});
  ar["checkpoint/operators"] << ops;
}

void QWL_SSE_Simulation::load(alps::hdf5::archive& ar) {
  alps::params p;
  uint64_t version,chain;
  std::vector<std::array<uint64_t,3>> graph;
  std::vector<double> factors;
  std::vector<uint32_t> ferro;
  ar["/parameters"] >> p;
  ar["checkpoint/version"] >> version;
  ar["checkpoint/chain"] >> chain;
  ar["checkpoint/bonds"] >> graph;
  ar["checkpoint/factors"] >> factors;
  ar["checkpoint/ferromagnetic"] >> ferro;
  if (checkpoint_parameters(p)!=checkpoint_parameters(parameters) || p.exists("SWEEPS")!=parameters.exists("SWEEPS") || version!=1 || chain!=chain_ || graph!=topology() || factors!=matrix_factor || ferro!=is_ferromagnetic)
    throw std::invalid_argument("QWL checkpoint model, topology or chain mismatch");
  auto saved=static_cast<qwl_state const&>(*this);
  ar["checkpoint/sweeps"] >> saved.sweeps;
  ar["checkpoint/production_sweeps"] >> saved.production_sweeps;
  ar["checkpoint/thistime"] >> saved.thistime;
  ar["checkpoint/non0"] >> saved.non0;
  ar["checkpoint/norder"] >> saved.norder;
  ar["checkpoint/logf_step"] >> saved.logf_step;
  ar["checkpoint/block_sweeps"] >> saved.block_sweeps;
  ar["checkpoint/block_sweeps_total"] >> saved.block_sweeps_total;
  ar["checkpoint/logf"] >> saved.logf;
  ar["checkpoint/minimum_histogram"] >> saved.minimum_histogram;
  ar["checkpoint/flatness_treshold"] >> saved.flatness_treshold;
  ar["checkpoint/doing_multicanonical"] >> saved.doing_multicanonical;
  ar["checkpoint/upwalker"] >> saved.upwalker;
  ar["checkpoint/all_done"] >> saved.all_done;
  ar["checkpoint/state"] >> saved.state;
  ar["checkpoint/g"] >> saved.g;
  ar["checkpoint/histo"] >> saved.histo;
  ar["checkpoint/histoup"] >> saved.histoup;
  ar["checkpoint/nmeasurements"] >> saved.nmeasurements;
  ar["checkpoint/uniform_structure_factor"] >> saved.uniform_structure_factor;
  ar["checkpoint/staggered_structure_factor"] >> saved.staggered_structure_factor;
  ar["checkpoint/transition_prob"] >> saved.transition_prob;
  std::vector<std::array<uint32_t,2>> ops;
  ar["checkpoint/operators"] >> ops;
  if (ops.size()!=L || saved.state.size()!=num_sites() || saved.non0!=saved.norder ||
      saved.production_sweeps>saved.sweeps || (!saved.doing_multicanonical && saved.production_sweeps) ||
      saved.norder<norder_min || saved.norder>norder_max || saved.sweeps>=UINT64_MAX/L ||
      !saved.logf_step || saved.logf_step>logf_steps_total || !saved.block_sweeps_total ||
      saved.block_sweeps>=saved.block_sweeps_total || saved.thistime>saved.sweeps*L ||
      saved.logf!=std::ldexp(initial_logf,1-int(saved.logf_step)) ||
      saved.minimum_histogram!=(use_zhou_bhatt ? std::ldexp(1/initial_logf,int(saved.logf_step)-1) : 0.) ||
      saved.flatness_treshold!=flatness_treshold || (saved.all_done && !saved.doing_multicanonical) ||
      (saved.doing_multicanonical && saved.logf_step!=logf_steps_total))
    throw std::invalid_argument("Invalid QWL checkpoint counters");
  auto propagated=saved.state;
  for (auto value:propagated) if (value>1) throw std::invalid_argument("Invalid QWL spin");
  unsigned nonidentity=0;
  for (size_t i=0;i<ops.size();++i) {
    auto const [op,b]=ops[i];
    if (op>3 || b>=num_bonds()) throw std::invalid_argument("Invalid QWL operator");
    saved.operator_string[i].op=op; saved.operator_string[i].bond=b;
    if (!op) continue;
    ++nonidentity;
    auto edge=bond(b);
    auto& left=propagated[source(edge)]; auto& right=propagated[target(edge)];
    if (op==1 && matrix_element(inhomogeneous_bond_type(edge),left,right)<=0)
      throw std::invalid_argument("Zero-weight QWL diagonal vertex");
    if (op>=2) {
      if (left!=(op==3) || right!=(op==2)) throw std::invalid_argument("Invalid QWL off-diagonal vertex");
      left^=1; right^=1;
    }
  }
  if (nonidentity!=saved.non0 || propagated!=saved.state)
    throw std::invalid_argument("QWL operator string is not periodic");
  for (auto const* histogram:{&saved.histo,&saved.histoup,&saved.nmeasurements,
                             &saved.uniform_structure_factor,&saved.staggered_structure_factor,&saved.transition_prob})
    if (histogram->size() && histogram->min()<0) throw std::invalid_argument("Negative QWL histogram");
  if (ar.list_children("measurements").size()!=measurements.size())
    throw std::invalid_argument("Unexpected QWL checkpoint measurements");
  alps::alea::hdf5_serializer codec(ar,"measurements");
  for (auto const& entry:measurements) {
    auto const& name=entry.first;
    std::visit([&](auto const& original) {
    using A=typename std::decay_t<decltype(original)>::element_type;
    A restored;
    alps::alea::deserialize(codec,ar.encode_segment(name),restored);
    if (restored.size()!=original->size() || restored.count()>saved.sweeps*L)
      throw std::invalid_argument("Invalid QWL measurement shape or count");
    if constexpr (std::is_same_v<A,alps::alea::batch_acc<double>>) {
      if (restored.num_batches()!=bins_ || restored.current_batch_size()!=restored.cursor().factor())
        throw std::invalid_argument("Invalid QWL timing batch layout");
    } else if constexpr (std::is_same_v<A,alps::alea::mean_acc<double>>) {
      uint64_t expected=saved.all_done;
      for (auto const& prefix:{std::string("Coefficients "),std::string("Total Sweeps ")})
        if (name.compare(0,prefix.size(),prefix)==0)
          expected=std::stoul(name.substr(prefix.size()))<saved.logf_step || saved.doing_multicanonical;
      if (restored.count()!=expected) throw std::invalid_argument("QWL phase/measurement mismatch");
    }
  },entry.second);
  }
  uint64_t expected_block=use_zhou_bhatt ? g.size() : uint64_t(parms.value_or_default("BLOCK_SWEEPS",10000));
  if (use_zhou_bhatt) for (unsigned step=1;step<saved.logf_step;++step)
    expected_block=static_cast<uint64_t>(expected_block*1.41);
  if (saved.doing_multicanonical && p.exists("SWEEPS")) expected_block=p["SWEEPS"].as<uint64_t>();
  if (saved.block_sweeps_total!=expected_block ||
      (saved.doing_multicanonical && saved.block_sweeps!=saved.production_sweeps%expected_block) ||
      (saved.doing_multicanonical && p.exists("SWEEPS") &&
       (saved.all_done ? saved.production_sweeps!=expected_block : saved.production_sweeps>=expected_block)))
    throw std::invalid_argument("QWL refinement/production schedule mismatch");
  bool extend=false;
  if (parameters.exists("SWEEPS") && saved.doing_multicanonical) {
    auto target=parameters["SWEEPS"].as<uint64_t>();
    if (saved.production_sweeps>target) throw std::invalid_argument("QWL production target precedes saved samples");
    extend=saved.all_done && saved.production_sweeps<target;
    saved.block_sweeps_total=target;
    if (extend) { saved.all_done=false; saved.block_sweeps=saved.production_sweeps; }
  }
  auto current=parameters;
  mcbase::load(ar);
  parameters=std::move(current);
  if (extend) for (auto const* name:{"Coefficients","Histogram","Fraction","Offset","Total Sweeps",
      "Uniform Structure Factor Coefficients","Staggered Structure Factor Coefficients"})
    if (measurements.count(name)) measurement<alps::alea::mean_acc<double>>(name)->reset();
  static_cast<qwl_state&>(*this)=std::move(saved);
  vertexbuild();
}

void QWL_SSE_Simulation::update() {
  if (all_done) throw std::logic_error("QWL run is complete");
  if (sweeps>=UINT64_MAX/L) throw std::overflow_error("QWL sweep counter overflow");
  ++sweeps;
  ++block_sweeps;
  if (doing_multicanonical)
    diagonal_update_multicanonical();
  else
    diagonal_update();
  if (non0) {
    vertexbuild();
    loopupdate();
  }
  else
    for (site_iterator it=sites().first;it!=sites().second;++it)
      state[*it]=random_state();
  if (doing_multicanonical) { ++production_sweeps; measure_magnetism(); }
  if (block_sweeps>=block_sweeps_total) {
    if (doing_multicanonical) {
      if (!all_done) {
        block_sweeps=0;
        if (histo.min()>=minimum_histogram || parms.defined("SWEEPS")) {
          store_histograms("multicanonical");
          cerr << "[" << norder_min << "-" << norder_max << "]  all done." << endl;
          all_done=1;
        }
      }
    }
    else {
      g.subtract();
      block_sweeps=0;
      double histomin=histo.min();
      double histoflatness=histo.flatness();
      cerr << "[" << norder_min << "-" << norder_max << "]  step "  << logf_step << ", ln[f]="<< logf
           << " : flatness="<< histoflatness << " ratio=" << histomin/minimum_histogram <<endl;
      if (histoflatness<flatness_treshold && histomin>=minimum_histogram) {
        store_histograms(std::to_string(logf_step));
        histo.fill(0.0);
        if (logf_step==logf_steps_total) {
          block_sweeps_total=parms.value_or_default("SWEEPS",block_sweeps_total);
          for (unsigned int i=g.left();i<g.right();++i)
            transition_prob[i]=exp(g[i]-g[i+1]);
          transition_prob[g.right()]=0.;
          cerr << "[" << norder_min << "-" << norder_max << "]  continuing using final weights..." << endl;
          doing_multicanonical=1;
        }
        else {
          ++logf_step;
          logf/=2.0;
          if (use_zhou_bhatt) {
            if (block_sweeps_total>UINT64_MAX/2) throw std::overflow_error("QWL refinement block overflow");
            block_sweeps_total=static_cast<uint64_t>(block_sweeps_total*1.41);
            minimum_histogram*=2;
          }
        }
      }
    }
  }
}


void QWL_SSE_Simulation::diagonal_update() { // fixed length (vector) representation in WL
  unsigned int Lmnon0=L-non0;
  for (vector<vertextype>::iterator operator_str_iterator=operator_string.begin();operator_str_iterator!=operator_string.end();++operator_str_iterator) {
    switch(operator_str_iterator->op) {
    // insert staebchenspiel here
    case 0 : // identity
      if (norder<g.right()) {
        operator_str_iterator->bond=random_int(0,num_bonds()-1);
        bond_descriptor this_bond=bond(operator_str_iterator->bond);
        double probability=exp(g[norder]-g[norder+1]);
        probability*=matrix_element(inhomogeneous_bond_type(this_bond),state[source(this_bond)],state[target(this_bond)]);
        probability*=num_bonds();
        probability*=beta;
        if (includeLfactor)
          probability/=Lmnon0;
        if (probability!=0. && ( (probability>=1) || (random_01()<probability) )) {
          operator_str_iterator->op=1;
          --Lmnon0;
          ++norder;
        }
      }
      break;
    case 1 : // diagonal
      if (norder>g.left()) {
        bond_descriptor this_bond=bond(operator_str_iterator->bond);
        double probability=exp(g[norder]-g[norder-1]);
        probability/=matrix_element(inhomogeneous_bond_type(this_bond),state[source(this_bond)],state[target(this_bond)]);
        probability/=num_bonds();
        probability/=beta;
        if (includeLfactor)
          probability*=Lmnon0+1;
        if (probability!=0. && ( (probability>=1) || (random_01()<probability) )) {
            operator_str_iterator->op=0;
          ++Lmnon0;
          --norder;
        }
      }
      break;
    default : // off-diagonal
      state[source(bond(operator_str_iterator->bond))]=get_state(operator_str_iterator->vs,2);
      state[target(bond(operator_str_iterator->bond))]=get_state(operator_str_iterator->vs,3);
    };
    if (norder==g.right()) upwalker=0;
    else if (norder==g.left())  upwalker=1;
    g[norder]+=logf;
    ++histo[norder];
  };
  non0=L-Lmnon0;
}


void QWL_SSE_Simulation::diagonal_update_multicanonical() {
  unsigned int Lmnon0=L-non0;
  for (vector<vertextype>::iterator operator_str_iterator=operator_string.begin();operator_str_iterator!=operator_string.end();++operator_str_iterator) {
    switch(operator_str_iterator->op) {
    case 0 :
      if (norder<g.right()) {
        operator_str_iterator->bond=random_int(0,num_bonds()-1);
        bond_descriptor this_bond=bond(operator_str_iterator->bond);
        double probability=transition_prob[norder];
        probability*=matrix_element(inhomogeneous_bond_type(this_bond),state[source(this_bond)],state[target(this_bond)]);
        probability*=num_bonds();
        probability*=beta;
        if (includeLfactor)
          probability/=Lmnon0;
        if (probability!=0. && ( (probability>=1) || (random_01()<probability) )) {
          operator_str_iterator->op=1;
          --Lmnon0;
          ++norder;
        }
      }
      break;
    case 1 :
      if (norder>g.left()) {
        bond_descriptor this_bond=bond(operator_str_iterator->bond);
        double probability=1./transition_prob[norder-1];
        probability/=matrix_element(inhomogeneous_bond_type(this_bond),state[source(this_bond)],state[target(this_bond)]);
        probability/=num_bonds();
        probability/=beta;
        if (includeLfactor)
          probability*=Lmnon0+1;
        if (probability!=0. && ( (probability>=1) || (random_01()<probability) )) {
            operator_str_iterator->op=0;
          ++Lmnon0;
          --norder;
        }
      }
      break;
    default :
      state[source(bond(operator_str_iterator->bond))]=get_state(operator_str_iterator->vs,2);
      state[target(bond(operator_str_iterator->bond))]=get_state(operator_str_iterator->vs,3);
    };
    ++thistime;
    if (upwalker && norder==g.right()) {
      record("Time Up",(double)thistime);
      record("Time Total",(double)thistime);
      thistime=0;
     upwalker=0;
    }
    else if (!upwalker && norder==g.left()) {
      record("Time Down",(double)thistime);
      record("Time Total",(double)thistime);
      thistime=0;
      upwalker=1;
    }
    ++histo[norder];
    if (upwalker) ++histoup[norder];
  };
  non0=L-Lmnon0;
}


void QWL_SSE_Simulation::vertexbuild() {
  for (unsigned int j=0;j<num_sites();++j)
    statev[j]=-(j+1);
  unsigned int i=0;
  for (vector <vertextype>::iterator operator_str_iterator=operator_string.begin();operator_str_iterator!=operator_string.end();++operator_str_iterator) {
    if (operator_str_iterator->op) {
      unsigned int bn0=source(bond(operator_str_iterator->bond));
      unsigned int bn1=target(bond(operator_str_iterator->bond));
      set_statevector_in(operator_str_iterator->vs,state[bn0],state[bn1]);
      operator_str_iterator->vv[0]=statev[bn0];
      operator_str_iterator->vv[1]=statev[bn1];
      operator_str_iterator->vn[0]=staten[bn0];
      operator_str_iterator->vn[1]=staten[bn1];
      if (statev[bn0] >= 0) {
        operator_string[statev[bn0]].vv[staten[bn0]]=i;
        operator_string[statev[bn0]].vn[staten[bn0]]=0;
      }
      if (statev[bn1] >= 0) {
        operator_string[statev[bn1]].vv[staten[bn1]]=i;
        operator_string[statev[bn1]].vn[staten[bn1]]=1;
      }
      if (operator_str_iterator->op!=1) {
        if (operator_str_iterator->op == 2) {
          ++state[bn0];
          --state[bn1];
        }
        else {
          --state[bn0];
          ++state[bn1];
        }
      }
      set_statevector_out(operator_str_iterator->vs,state[bn0],state[bn1]);
      statev[bn0]=i;
      statev[bn1]=i;
      staten[bn0]=2;
      staten[bn1]=3;
    }
    ++i;
  };
  i=0;
  for (vector <vertextype>::iterator operator_str_iterator=operator_string.begin();operator_str_iterator!=operator_string.end();++operator_str_iterator) {
    if (operator_str_iterator->op) {
      int j=-(operator_str_iterator->vv[0]+1);
      if (j >= 0) {
        operator_str_iterator->vv[0]=statev[j];
        operator_str_iterator->vn[0]=staten[j];
        operator_string[statev[j]].vv[staten[j]]=i;
        operator_string[statev[j]].vn[staten[j]]=0;
      }
      j=-(operator_string[i].vv[1]+1);
      if (j >= 0) {
        operator_str_iterator->vv[1]=statev[j];
        operator_str_iterator->vn[1]=staten[j];
        operator_string[statev[j]].vv[staten[j]]=i;
        operator_string[statev[j]].vn[staten[j]]=1;
      }
      for (int l=0;l<4;++l)
        operator_str_iterator->vvisited[l]=0;
    }
    else
      for (int l=0;l<4;++l)
        operator_str_iterator->vvisited[l]=1;
    ++i;
  };
}

void QWL_SSE_Simulation::loopupdate() {
  unsigned int i0;
  unsigned int ir;
  boost::uint32_t n0;
  boost::uint32_t nr;
  boost::uint32_t nex;
  i0=0;
  n0=0;
  int flipflag;
  do {
    while ( ( (i0!=L) ?  operator_string[i0].vvisited[n0] : 0) ) {
      if (n0==3) {
        ++i0;
        n0=0;
      }
      else
       ++n0;
    };
    if (i0==L) break;
    flipflag=(random_01()>0.5);
    ir=i0;
    nr=n0;
    do {
      operator_string[ir].vvisited[nr]=1;
      if (flipflag)
        transform(operator_string[ir].vs,nr);
      nex=exit_leg(operator_string[ir],nr);
      operator_string[ir].vvisited[nex]=1;
      if (flipflag)
        transform(operator_string[ir].vs,nex);
      if (ir==i0 && nex==n0)
        break;
      nr=operator_string[ir].vn[nex];
      ir=operator_string[ir].vv[nex];
    } while((ir!=i0)||(nr!=n0));
  } while (1);
  for (vector <vertextype>::iterator operator_str_iterator=operator_string.begin();operator_str_iterator!=operator_string.end();++operator_str_iterator)
    if (operator_str_iterator->op)
      switch (get_state(operator_str_iterator->vs,2)-get_state(operator_str_iterator->vs,0)) {
        case 0:
          operator_str_iterator->op=1;
          break;
        case 1:
          operator_str_iterator->op=2;
          break;
        default:
          operator_str_iterator->op=3;
      }
  for (unsigned int j=0;j<num_sites();++j)
    state[j]= (statev[j]>=0) ? get_state(operator_string[statev[j]].vs,staten[j]) : random_state();
}


inline boost::uint32_t QWL_SSE_Simulation::exit_leg(vertextype&  opn, boost::uint32_t nr) {
  if (is_ferromagnetic[inhomogeneous_bond_type(bond(opn.bond))])
    return 3-nr;
  if (!nr)
    return 1;
  if (nr==1)
    return 0;
  if (nr==2)
    return 3;
  return 2;
}


inline double QWL_SSE_Simulation::matrix_element(const int b_t, const int s0, const int s1) {
  if (is_ferromagnetic[b_t])
    return ( (s0==s1) ? matrix_factor[b_t] : 0. );
  return ( (s0==s1) ? 0. : matrix_factor[b_t] );
}


void QWL_SSE_Simulation::init_tables() {
  for (auto [it,end]=sites();it!=end;++it) {
    auto type=inhomogeneous_site_type(*it);
    auto local=alps::get_matrix(double(),model().site_term(type),model().basis().site_basis(type),parms);
    for (auto value:std::vector<double>(local.data(),local.data()+local.num_elements()))
      if (!std::isfinite(value) || std::abs(value)>1e-12)
        throw std::invalid_argument("QWL supports spin-1/2 isotropic exchange without on-site terms");
  }
  int max_bond_type=0;
  for (bond_iterator it=bonds().first; it!=bonds().second;++it) {
    int this_bond_type=inhomogeneous_bond_type(*it);
    if (this_bond_type>=max_bond_type) max_bond_type=this_bond_type;
  }
  is_ferromagnetic.resize(max_bond_type+1);
  matrix_factor.resize(max_bond_type+1);
  offset=0;
  for (bond_iterator it=bonds().first; it!=bonds().second;++it) {
    if (source(*it)==target(*it)) throw std::invalid_argument("QWL requires distinct bond endpoints");
    Parameters p(parms);
    if (inhomogeneous())
      throw_if_xyz_defined(parms,*it);
    if (inhomogeneous_sites()) {
      p << coordinate_as_parameter(source(*it));
      p << coordinate_as_parameter(target(*it));
    }
    if (inhomogeneous_bonds())
      p << coordinate_as_parameter(*it);
    int this_bond_type=inhomogeneous_bond_type(*it);
    int source_site_type=inhomogeneous_site_type(source(*it));
    int target_site_type=inhomogeneous_site_type(target(*it));
    int source_num_states = model().basis().site_basis(source_site_type).num_states();
    int target_num_states = model().basis().site_basis(target_site_type).num_states();
    if (source_num_states!=2 || target_num_states!=2)
      boost::throw_exception(std::runtime_error("This model cannot be simulated with this code"));
    boost::multi_array<double,2> source_Szmatrix_symbolic =
      alps::get_matrix(double(),
      SiteOperator("Sz"),
      model().basis().site_basis(source_site_type),
      p);
    boost::multi_array<double,2> target_Szmatrix_symbolic =
      alps::get_matrix(double(),
      SiteOperator("Sz"),
      model().basis().site_basis(target_site_type),
      p);
    int source_spin_up_state=(source_Szmatrix_symbolic[0][0]>0.) ? 0 : 1;
    int target_spin_up_state=(target_Szmatrix_symbolic[0][0]>0.) ? 0 : 1;
    boost::multi_array<double,4> bondhamiltonian =
      alps::get_matrix(double(),
      model().bond_term(this_bond_type),
      model().basis().site_basis(source_site_type),
      model().basis().site_basis(target_site_type),
      p);
    double parallel=bondhamiltonian[source_spin_up_state][target_spin_up_state][source_spin_up_state][target_spin_up_state];
    for (int a=0;a<2;++a) for (int b=0;b<2;++b)
      for (int c=0;c<2;++c) for (int d=0;d<2;++d) {
        double expected=(a==c && b==d) ? (a==b ? parallel : -parallel) :
                        (a!=b && a==d && b==c ? 2*parallel : 0.);
        double actual=bondhamiltonian[a][b][c][d];
        if (!std::isfinite(actual) || std::abs(actual-expected)>1e-12*std::max(1.,std::abs(parallel)))
          throw std::invalid_argument("QWL requires isotropic spin-1/2 exchange");
      }
    is_ferromagnetic[this_bond_type]=
      bondhamiltonian[source_spin_up_state][target_spin_up_state]
                     [source_spin_up_state][target_spin_up_state]<0? 1 : 0;
    matrix_factor[this_bond_type]=
      fabs(bondhamiltonian[source_spin_up_state][target_spin_up_state]
                          [source_spin_up_state][target_spin_up_state])*2.;
    offset+=
      fabs(bondhamiltonian[source_spin_up_state][target_spin_up_state]
                          [source_spin_up_state][target_spin_up_state]);
  }
  if (alps::has_sign_problem(model(),static_cast<graph_helper<> const&>(*this),parms))
    throw std::invalid_argument("QWL cannot sample a model with a sign problem");
  if (std::none_of(matrix_factor.begin(),matrix_factor.end(),[](double x){ return x>0; }))
    throw std::invalid_argument("QWL requires a nonzero exchange coupling");
}


void QWL_SSE_Simulation::init_measurements() {
  auto add=[&](std::string name,size_t size) {
    measurements.emplace(std::move(name),std::make_shared<alps::alea::mean_acc<double>>(size));
  };
  add("Offset",1);
  const auto first=std::min<unsigned>(parms.value_or_default("START_STORING",logf_steps_total),logf_steps_total);
  for (unsigned pos=first;pos<=logf_steps_total;++pos) {
    add("Coefficients "+std::to_string(pos),g.size());
    add("Total Sweeps "+std::to_string(pos),1);
  }
  for (auto name:{"Coefficients","Histogram","Fraction"}) add(name,g.size());
  add("Total Sweeps",1);
  for (auto name:{"Time Up","Time Down","Time Total"}) {
    measurements.emplace(name,std::make_shared<alps::alea::batch_acc<double>>(1,bins_));
    measurements.emplace("Autocorrelation: "+std::string(name),std::make_shared<alps::alea::autocorr_acc<double>>(1));
  }
  if (measure_magnetics) {
    add("Uniform Structure Factor Coefficients",g.size());
    if (is_bipartite()) add("Staggered Structure Factor Coefficients",g.size());
  }
}


void QWL_SSE_Simulation::measure_magnetism() {
  if (measure_magnetics) {
    ++nmeasurements[norder];
    double mag=0;
    double smag=0;
    for (unsigned int i=0;i<num_sites();++i) {
      double local_state=decode_spin(state[i]);
      mag+=local_state;
      if (is_bipartite())
        smag+=parity(i)*local_state;
    }
    uniform_structure_factor[norder]+=sqr(mag);
    if (is_bipartite())
      staggered_structure_factor[norder]+=sqr(smag);
  }
}


void QWL_SSE_Simulation::store_histograms(string pos) {
  if (pos=="multicanonical") {
    valarray<double> gval=g.getvalarray(norder_min,norder_max);
    valarray<double> histoval=histo.getvalarray(norder_min,norder_max);
    valarray<double> histoupval=histoup.getvalarray(norder_min,norder_max);
    histoupval/=histoval;
    histoval/=block_sweeps_total;
    double histovalsum=histoval.sum();
    histoval/=histovalsum;
    if (histoval[0]<=0) throw std::runtime_error("QWL production did not visit its normalization order; increase SWEEPS");
    double red=gval[0]+log(histoval[0])-(norder_min==0 ? num_sites()*log(2.) : 0.);
    for (int i=0;i<gval.size();++i)
      gval[i]+=log(histoval[i])-red;
    if (!includeLfactor) for (unsigned i=0;i<gval.size();++i)
      gval[i]+=std::lgamma(double(L-norder_min-i)+1)-std::lgamma(double(L-norder_min)+1);
    record("Coefficients",gval);
    record("Histogram",histoval);
    record("Fraction",histoupval);
    record("Total Sweeps",(double)sweeps);
    record("Offset",offset);
    if (measure_magnetics) {
      valarray<double> nmeasurementsval=nmeasurements.getvalarray(norder_min,norder_max);
      valarray<double> uniform_structure_factor_val=uniform_structure_factor.getvalarray(norder_min,norder_max);
      uniform_structure_factor_val/=nmeasurementsval;
      record("Uniform Structure Factor Coefficients",uniform_structure_factor_val);
      if (is_bipartite()) {
        valarray<double> staggered_structure_factor_val=staggered_structure_factor.getvalarray(norder_min,norder_max);
        staggered_structure_factor_val/=nmeasurementsval;
        record("Staggered Structure Factor Coefficients",staggered_structure_factor_val);
      }
    }
  }
  else {
    int posmin=parms.value_or_default("START_STORING",logf_steps_total);
    if (posmin>logf_steps_total)
      posmin=logf_steps_total;
    if (atoi(pos.c_str())>=posmin) {
      valarray<double> gval=g.getvalarray(norder_min,norder_max);
      double red=gval[0]-(norder_min==0 ? num_sites()*log(2.) : 0.);
      gval-=red;
      if (!includeLfactor) for (unsigned i=0;i<gval.size();++i)
        gval[i]+=std::lgamma(double(L-norder_min-i)+1)-std::lgamma(double(L-norder_min)+1);
      record("Coefficients "+pos,gval);
      record("Total Sweeps "+pos,(double)sweeps);
    }
  }
}

#endif
