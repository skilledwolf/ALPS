/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1997-2010 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#pragma once
#include <type_traits>
#include "loop_config.h"
#include <looper/cluster.h>
#include <looper/operator.h>
#include <looper/permutation.h>
#include <looper/temperature.h>
#include <looper/type.h>
#include <looper/union_find.h>
#include "../simulation.hpp"
#include <alps/hdf5/stdarray.hpp>


namespace looper {

template<class MC>
class loop_worker : public native_qmc::simulation, private loop_config {
public:
  typedef MC mc_type;
  static constexpr bool continuous_time = std::is_same_v<MC, looper::path_integral>;

  using time_t = std::conditional_t<continuous_time, loop_config::time_t, int>;
  typedef looper::local_operator<mc_type, loop_graph_t, time_t> local_operator_t;
  typedef std::vector<local_operator_t> operator_string_t;
  typedef typename operator_string_t::iterator operator_iterator;

  typedef looper::union_find::node cluster_fragment_t;
  typedef looper::cluster_info cluster_info_t;

  typedef typename looper::estimator<measurement_set, mc_type, lattice_t, time_t>::type estimator_t;
  using weight_parameter_type=Eigen::Vector2d;

  loop_worker(alps::params const& p,size_t bins,size_t chain,double initial_beta=0.);
  void init_observables(native_qmc::simulation& obs);

  bool is_thermalized() const { return steps_>=parameters["THERMALIZATION"].template as<uint64_t>(); }
  double fraction_completed() const override {
    return is_thermalized() ? double(steps_-parameters["THERMALIZATION"].template as<uint64_t>())/parameters["SWEEPS"].template as<uint64_t>() : 0.;
  }
  uint64_t completed_sweeps() const { return steps_; }
  size_t site_count() const { return num_sites(lattice.rg()); }
  double volume() const { return lattice.volume(); }
  void update() override {
    record_measurements(is_thermalized());
    run(*this,1./temperature(steps_+1));
  }
  void measure() override {}
  void run(native_qmc::simulation& obs,double inverse_temperature);

  // Called after a completed sweep; field_energy_ is reconstructed each update.
  weight_parameter_type weight_parameter() const { return {double(operators.size()),field_energy_}; }
  static double log_weight(weight_parameter_type const& weight,double inverse_temperature) {
    return std::log(inverse_temperature)*weight[0]-inverse_temperature*weight[1];
  }
  void save(alps::hdf5::archive& ar) const override;
  void load(alps::hdf5::archive& ar) override;


protected:
  std::vector<double> model_state() const {
    std::vector<double> result{model.graph_weight(),model.energy_offset(),lattice.volume()};
    for (double value:model.graph_weights()) result.push_back(value);
    for (double value:model.field()) result.push_back(value);
    for (int value:model.site_sign()) result.push_back(value);
    for (int value:model.bond_sign()) result.push_back(value);
    for (auto [it,end]=sites(lattice.rg());it!=end;++it) {
      auto [first,last]=sites(lattice,*it); result.push_back(std::distance(first,last));
    }
    for (auto [it,end]=bonds(lattice.vg());it!=end;++it) {
      result.push_back(source(*it,lattice.vg())); result.push_back(target(*it,lattice.vg()));
    }
    return result;
  }
  void build();
  void connect(local_operator_t&);

  template<typename FIELD, typename SIGN, typename IMPROVE>
  void flip(native_qmc::simulation& obs);

private:
  // helpers
  lattice_t lattice;
  model_t model;

  // parameters
  looper::temperature temperature;
  double beta;
  bool use_improved_estimator;

  // configuration (checkpoint)
  uint64_t steps_=0;
  size_t chain_;
  std::vector<int> spins;
  std::vector<local_operator_t> operators;

  // observables
  double sign;
  double field_energy_=0;
  estimator_t estimator;

  // working vectors
  std::vector<int> spins_c;
  std::vector<local_operator_t> operators_p;
  std::vector<cluster_fragment_t> fragments;
  std::vector<int> current;
  std::vector<bool> to_flip;
  std::vector<cluster_info_t> clusters;
  std::vector<typename looper::estimate<estimator_t>::type> estimates;
  std::vector<int> perm;
};


//
// member functions of loop_worker
//

template<class MC>
loop_worker<MC>::loop_worker(alps::params const& p,size_t bins,size_t chain,double initial_beta)
  : native_qmc::simulation(p,bins,chain), lattice(alps::Parameters(alps::seed_disorder(p))),
    model(alps::Parameters(p), lattice, continuous_time),
    temperature([&] {
      alps::Parameters values(p);
      if (initial_beta>0) values["T"]=1/initial_beta;
      return values;
    }()), chain_(chain) {

  if (initial_beta>0 && temperature.annealing_steps())
    throw std::invalid_argument("Annealing and replica exchange cannot be combined");

  if (!site_count() || bins<2 || bins%2 || !std::isfinite(temperature.final()) ||
      temperature.final()<=0 || !std::isfinite(temperature.initial()) || temperature.initial()<=0)
    throw std::invalid_argument("Loop requires a nonempty lattice, positive finite temperature and even batch capacity >= 2");
  if (temperature.annealing_steps() > p["THERMALIZATION"].as<uint64_t>())
    boost::throw_exception(std::invalid_argument("longer annealing steps than thermalization"));

  model.check_parameter(continuous_time && support_longitudinal_field, support_negative_sign);

  use_improved_estimator = (!model.has_field()) && (!p.value_or("DISABLE_IMPROVED_ESTIMATOR",false));
  if (!use_improved_estimator) std::cout << "WARNING: improved estimator is disabled\n";

  // configuration
  int nvs = num_sites(lattice.vg());
  spins.resize(nvs); std::fill(spins.begin(), spins.end(), 0 /* all up */);
  spins_c.resize(nvs);
  current.resize(nvs);
  perm.resize(max_virtual_sites(lattice));

  // initialize estimators
  estimator.initialize(alps::Parameters(p), lattice, model.is_signed(), use_improved_estimator);
  is_signed_=model.is_signed();
  init_observables(*this);
}

template<class MC>
void loop_worker<MC>::init_observables(native_qmc::simulation& obs) {
  obs.add_measurement("Temperature",1,false);
  obs.add_measurement("Inverse Temperature",1,false);
  obs.add_measurement("Volume",1,false);
  obs.add_measurement("Number of Sites",1,false);
  obs.add_measurement("Number of Clusters",1,false);
  if (model.is_signed()) {
    obs.add_measurement("Sign",1,false);
    if (use_improved_estimator) {
      obs.add_measurement("Weight of Zero-Meron Sector",1,false);
      obs.add_measurement("Sign in Zero-Meron Sector",1,false);
    }
  }
  looper::energy_estimator::init_observables(obs, model.is_signed());
  estimator.init_observables(obs, model.is_signed());
}

template<class MC>
void loop_worker<MC>::run(native_qmc::simulation& obs,double inverse_temperature) {
  ++steps_;
  beta = inverse_temperature;

  build();

  //   FIELD               SIGN                IMPROVE
  flip<boost::mpl::true_,  boost::mpl::true_,  boost::mpl::true_ >(obs);
  flip<boost::mpl::true_,  boost::mpl::true_,  boost::mpl::false_>(obs);
  flip<boost::mpl::true_,  boost::mpl::false_, boost::mpl::true_ >(obs);
  flip<boost::mpl::true_,  boost::mpl::false_, boost::mpl::false_>(obs);
  flip<boost::mpl::false_, boost::mpl::true_,  boost::mpl::true_ >(obs);
  flip<boost::mpl::false_, boost::mpl::true_,  boost::mpl::false_>(obs);
  flip<boost::mpl::false_, boost::mpl::false_, boost::mpl::true_ >(obs);
  flip<boost::mpl::false_, boost::mpl::false_, boost::mpl::false_>(obs);
}


//
// diagonal update and cluster construction
//

template<class MC>
void loop_worker<MC>::build() {
  // initialize spin & operator information
  std::copy(spins.begin(), spins.end(), spins_c.begin());
  std::swap(operators, operators_p); operators.resize(0);

  // initialize cluster information (setup cluster fragments)
  int nvs = num_sites(lattice.vg());
  fragments.resize(0); fragments.resize(nvs);
  for (int s = 0; s < nvs; ++s) current[s] = s;

  if constexpr (continuous_time) {
    auto r_time=[&] {
      double rate=beta*model.graph_weight();
      return rate>0 ? -std::log1p(-random())/rate : std::numeric_limits<double>::infinity();
    };
    double t = r_time();
    for (operator_iterator opi = operators_p.begin(); t < 1 || opi != operators_p.end();) {

      // diagonal update & labeling
      if (opi == operators_p.end() || t < opi->time()) {
        loop_graph_t g = model.choose_graph(random);
        if ((is_bond(g) && is_compatible(g, spins_c[source(pos(g), lattice.vg())],
                                            spins_c[target(pos(g), lattice.vg())])) ||
            (is_site(g) && is_compatible(g, spins_c[pos(g)]))) {
          operators.push_back(local_operator_t(g, t));
          t += r_time();
        } else {
          t += r_time();
          continue;
        }
      } else {
        if (opi->is_diagonal()) {
          ++opi;
          continue;
        } else {
          operators.push_back(*opi);
          ++opi;
        }
      }

      connect(operators.back());
    }
  } else {
    int nop = operators_p.size();
    double bw = beta * model.graph_weight();
    bool try_gap = true;
    for (operator_iterator opi = operators_p.begin(); try_gap || opi != operators_p.end();) {

      // diagonal update & labeling
      if (try_gap) {
        if ((nop+1) * random() < bw) {
          loop_graph_t g = model.choose_graph(random);
          if ((is_bond(g) && is_compatible(g, spins_c[source(pos(g), lattice.vg())],
                                              spins_c[target(pos(g), lattice.vg())])) ||
              (is_site(g) && is_compatible(g, spins_c[pos(g)]))) {
            operators.push_back(local_operator_t(g));
            ++nop;
          } else {
            try_gap = false;
            continue;
          }
        } else {
          try_gap = false;
          continue;
        }
      } else {
        if (opi->is_diagonal()) {
          if (bw * random() < nop) {
            --nop;
            ++opi;
            continue;
          } else {
            if (opi->is_site()) {
              opi->assign_graph(model.choose_diagonal(random, opi->loc(),
                spins_c[opi->pos()]));
            } else {
              opi->assign_graph(model.choose_diagonal(random, opi->loc(),
                spins_c[source(opi->pos(), lattice.vg())],
                spins_c[target(opi->pos(), lattice.vg())]));
            }
          }
        } else {
          if (opi->is_bond())
            opi->assign_graph(model.choose_offdiagonal(random, opi->loc(),
              spins_c[source(opi->pos(), lattice.vg())],
              spins_c[target(opi->pos(), lattice.vg())]));
        }
        operators.push_back(*opi);
        ++opi;
        try_gap = true;
      }

      connect(operators.back());
    }
  }

  // symmetrize spins
  if (max_virtual_sites(lattice) == 1) {
    for (int i = 0; i < nvs; ++i) unify(fragments, i, current[i]);
  } else {
    BOOST_FOREACH(looper::real_site_descriptor<lattice_t>::type rs, sites(lattice.rg())) {
      looper::virtual_site_iterator<lattice_t>::type vsi, vsi_end;
      boost::tie(vsi, vsi_end) = sites(lattice, rs);
      int offset = *vsi;
      int s2 = *vsi_end - *vsi;
      for (int i = 0; i < s2; ++i) perm[i] = i;
      looper::partitioned_random_shuffle(perm.begin(), perm.begin() + s2,
        spins.begin() + offset, spins_c.begin() + offset, random);
      for (int i = 0; i < s2; ++i) unify(fragments, offset+i, current[offset+perm[i]]);
    }
  }
}


template<class MC>
void loop_worker<MC>::connect(local_operator_t& op) {
  if (op.is_bond()) {
    int s0 = source(op.pos(), lattice.vg());
    int s1 = target(op.pos(), lattice.vg());
    if (op.is_offdiagonal()) {
      if constexpr (continuous_time)
        op.assign_graph(model.choose_offdiagonal(random, op.loc(),
        spins_c[s0], spins_c[s1]));
      spins_c[s0] ^= 1;
      spins_c[s1] ^= 1;
    }
    boost::tie(current[s0], current[s1], op.loop0, op.loop1) =
      reconnect(fragments, op.graph(), current[s0], current[s1]);
  } else {
    int s = op.pos();
    if (op.is_offdiagonal()) spins_c[s] ^= 1;
    boost::tie(current[s], op.loop0, op.loop1) = reconnect(fragments, op.graph(), current[s]);
  }
}

//
// cluster flip
//

template<class MC>
template<typename FIELD, typename SIGN, typename IMPROVE>
void loop_worker<MC>::flip(native_qmc::simulation& obs) {
  if (model.has_field() != FIELD() ||
      model.is_signed() != SIGN() ||
      use_improved_estimator != IMPROVE()) return;

  int nvs = num_sites(lattice.vg());

  // assign cluster id
  int nc = 0;
  BOOST_FOREACH(cluster_fragment_t& f, fragments) if (f.is_root()) f.set_id(nc++);
  BOOST_FOREACH(cluster_fragment_t& f, fragments) f.set_id(cluster_id(fragments, f));
  to_flip.resize(nc);
  clusters.resize(0); clusters.resize(nc);

  std::copy(spins.begin(), spins.end(), spins_c.begin());
  cluster_info_t::accumulator<cluster_fragment_t, FIELD, SIGN, IMPROVE>
    weight(clusters, fragments, model.field(), model.bond_sign(), model.site_sign());
  looper::accumulator<estimator_t, cluster_fragment_t, IMPROVE>
    accum(estimates, nc, lattice, estimator, fragments);
  for (unsigned int s = 0; s < nvs; ++s) {
    weight.start_bottom(s, time_t(0), s, spins_c[s]);
    accum.start_bottom(s, time_t(0), s, spins_c[s]);
  }
  time_t t = 0;
  int negop = 0; // number of operators with negative weights
  BOOST_FOREACH(local_operator_t& op, operators) {
    if constexpr (continuous_time) t=op.time();
    if (op.is_bond()) {
      if (!op.is_frozen_bond_graph()) {
        int b = op.pos();
        int s0 = source(b, lattice.vg());
        int s1 = target(b, lattice.vg());
        weight.end_b(op.loop_l0(), op.loop_l1(), t, b, s0, s1, spins_c[s0], spins_c[s1]);
        accum.end_b(op.loop_l0(), op.loop_l1(), t, b, s0, s1, spins_c[s0], spins_c[s1]);
        if (op.is_offdiagonal()) {
          spins_c[s0] ^= 1;
          spins_c[s1] ^= 1;
          if (SIGN()) negop += model.bond_sign(op.pos());
        }
        weight.begin_b(op.loop_u0(), op.loop_u1(), t, b, s0, s1, spins_c[s0], spins_c[s1]);
        accum.begin_b(op.loop_u0(), op.loop_u1(), t, b, s0, s1, spins_c[s0], spins_c[s1]);
      }
    } else {
      if (!op.is_frozen_site_graph()) {
        int s = op.pos();
        weight.end_s(op.loop_l(), t, s, spins_c[s]);
        accum.end_s(op.loop_l(), t, s, spins_c[s]);
        if (op.is_offdiagonal()) {
          spins_c[s] ^= 1;
          if (SIGN()) negop += model.site_sign(op.pos());
        }
        weight.begin_s(op.loop_u(), t, s, spins_c[s]);
        accum.begin_s(op.loop_u(), t, s, spins_c[s]);
      }
    }
    if constexpr (!continuous_time) ++t;
  }
  for (unsigned int s = 0; s < nvs; ++s) {
    weight.stop_top(current[s], time_t(continuous_time ? 1 : operators.size()), s, spins_c[s]);
    accum.stop_top(current[s], time_t(continuous_time ? 1 : operators.size()), s, spins_c[s]);
  }
  sign = ((negop & 1) == 1) ? -1 : 1;

  // accumulate cluster properties
  typename looper::collector<estimator_t>::type coll = get_collector(estimator);
  coll.set_num_operators(operators.size());
  coll.set_num_clusters(nc);
  if (IMPROVE()) {
    BOOST_FOREACH(typename looper::estimate<estimator_t>::type const& est, estimates) { coll += est; }
  }

  // determine whether clusters are flipped or not
  double improved_sign = sign;
  for (unsigned int c = 0; c < clusters.size(); ++c) {
    to_flip[c] = ((2*random()-1) < (FIELD() ? std::tanh(beta * clusters[c].weight) : 0));
    if (SIGN() && IMPROVE() && (clusters[c].sign & 1) == 1) improved_sign = 0;
  }

  // improved measurement
  if (IMPROVE())
    estimator.improved_measurement(obs, lattice, beta, improved_sign, spins, operators,
      spins_c, fragments, coll);

  // Normal estimators below observe the flipped configuration, so carry its
  // sign too. Improved estimators above retain their pre-flip cluster sign.
  BOOST_FOREACH(local_operator_t& op, operators)
    if (to_flip[fragments[op.loop_0()].id()] ^ to_flip[fragments[op.loop_1()].id()]) {
      op.flip();
      if (SIGN() && (op.is_bond() ? model.bond_sign(op.pos()) : model.site_sign(op.pos()))) sign=-sign;
    }
  for (int s = 0; s < nvs; ++s) if (to_flip[fragments[s].id()]) spins[s] ^= 1;

  //
  // measurement
  //

  obs.record("Temperature", 1/beta, 1.);
  obs.record("Inverse Temperature", beta, 1.);
  obs.record("Volume", (double)lattice.volume(), 1.);
  obs.record("Number of Sites", (double)num_sites(lattice.rg()), 1.);
  obs.record("Number of Clusters", coll.num_clusters(), 1.);

  // sign
  if (SIGN()) {
    if (IMPROVE()) {
      obs.record("Sign", improved_sign, 1.);
      if (alps::numeric::is_zero(improved_sign)) {
        obs.record("Weight of Zero-Meron Sector", 0., 1.);
      } else {
        obs.record("Weight of Zero-Meron Sector", 1., 1.);
        obs.record("Sign in Zero-Meron Sector", improved_sign, 1.);
      }
    } else {
      obs.record("Sign", sign, 1.);
    }
  }

  // energy
  double nop = coll.num_operators();
  field_energy_=0;
  double ene = model.energy_offset() - nop / beta;
  if (FIELD())
    for (unsigned int c = 0; c < clusters.size(); ++c)
      field_energy_ += (to_flip[c] ? -clusters[c].weight : clusters[c].weight);
  ene+=field_energy_;
  looper::energy_estimator::measurement(obs, lattice, beta, nop, sign, ene);

  // normal measurement
  estimator.normal_measurement(obs, lattice, beta, sign, spins, operators, spins_c);
}

template<class MC>
void loop_worker<MC>::save(alps::hdf5::archive& ar) const {
  alps::mcbase::save(ar);
  ar["checkpoint/version"] << uint64_t(1);
  ar["checkpoint/steps"] << steps_;
  ar["checkpoint/chain"] << uint64_t(chain_);
  ar["checkpoint/model"] << model_state();
  ar["checkpoint/spins"] << spins;
  std::vector<std::array<int,3>> ops;
  std::vector<double> times;
  for (auto const& op:operators) {
    ops.push_back({op.type()|(op.graph_type()<<2),op.is_bond()?1:0,op.pos()});
    if constexpr (continuous_time) times.push_back(op.time());
  }
  ar["checkpoint/operators"] << ops;
  if constexpr (continuous_time) ar["checkpoint/times"] << times;
}

template<class MC>
void loop_worker<MC>::load(alps::hdf5::archive& ar) {
  uint64_t version,steps,chain;
  std::vector<double> fingerprint,times;
  std::vector<int> state;
  std::vector<std::array<int,3>> ops;
  ar["checkpoint/version"] >> version;
  ar["checkpoint/steps"] >> steps;
  ar["checkpoint/chain"] >> chain;
  ar["checkpoint/model"] >> fingerprint;
  ar["checkpoint/spins"] >> state;
  ar["checkpoint/operators"] >> ops;
  if constexpr (continuous_time) ar["checkpoint/times"] >> times;
  if (version!=1 || chain!=chain_ || fingerprint!=model_state() || state.size()!=spins.size() ||
      !std::all_of(state.begin(),state.end(),[](int spin){return spin==0 || spin==1;}) ||
      (continuous_time && times.size()!=ops.size()))
    throw std::invalid_argument("Invalid loop checkpoint model or shape");
  operator_string_t restored;
  auto propagated=state;
  for (size_t i=0;i<ops.size();++i) {
    auto [type,bond,pos]=ops[i];
    int graph=type>>2;
    if (type<0 || (type&3)>1 || (bond!=0 && bond!=1) || pos<0 ||
        (bond && (pos>=num_bonds(lattice.vg()) || !loop_graph_t::bond_graph_t::is_valid_gid(graph))) ||
        (!bond && (pos>=state.size() || !loop_graph_t::site_graph_t::is_valid_gid(graph))))
      throw std::invalid_argument("Invalid loop checkpoint operator");
    auto location=bond ? loop_graph_t::location_t::bond_location(pos) : loop_graph_t::location_t::site_location(pos);
    if constexpr (continuous_time) {
      if (!std::isfinite(times[i]) || times[i]<0 || times[i]>=1 || (i && times[i]<=times[i-1]))
        throw std::invalid_argument("Invalid loop checkpoint time ordering");
      restored.emplace_back(type,location,times[i]);
    } else restored.emplace_back(type,location);
    auto const g=restored.back().graph();
    // Cluster flips retain the auxiliary graph label until the next build.
    // Its diagonal compatibility predicate need not hold after an operator
    // becomes offdiagonal; the XXZ matrix element exchanges antiparallel spins.
    const bool compatible=bond && (type&1)
      ? !restored.back().is_frozen_bond_graph() && propagated[source(pos,lattice.vg())]!=propagated[target(pos,lattice.vg())]
      : bond ? is_compatible(g,propagated[source(pos,lattice.vg())],propagated[target(pos,lattice.vg())])
             : is_compatible(g,propagated[pos]);
    if (!compatible)
      throw std::invalid_argument("Loop checkpoint operator incompatible with worldline spins");
    if (type&1) {
      if (bond) { propagated[source(pos,lattice.vg())]^=1; propagated[target(pos,lattice.vg())]^=1; }
      else propagated[pos]^=1;
    }
  }
  for (auto [it,end]=sites(lattice.rg());it!=end;++it) {
    auto [first,last]=sites(lattice,*it);
    int difference=0;
    for (;first!=last;++first) difference+=state[*first]-propagated[*first];
    if (difference) throw std::invalid_argument("Nonperiodic loop checkpoint worldlines");
  }
  validate_measurements(ar,steps);
  auto requested=parameters;
  alps::mcbase::load(ar);
  parameters=std::move(requested); // retain the requested extended SWEEPS limit
  steps_=steps; spins=std::move(state); operators=std::move(restored);
}

} // namespace looper
