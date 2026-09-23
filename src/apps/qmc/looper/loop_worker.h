#ifndef ALPS_LOOP_WORKER_H
#define ALPS_LOOP_WORKER_H

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

#include "loop_config.h"
#include <looper/cluster.h>
#include <looper/evaluator_impl.h>
#include <looper/montecarlo.h>
#include <looper/operator.h>
#include <looper/permutation.h>
#include <looper/temperature.h>
#include <looper/type.h>
#include <looper/union_find.h>
#ifdef HAVE_PARAPACK_13
# include <alps/parapack/serial.h>
#else
# include <alps/parapack/worker.h>
#endif
#include <alps/parapack/exchange.h>

namespace looper {

// Representation differences used by the shared cluster and measurement workflow.
template<class Representation> struct worker_representation;

template<> struct worker_representation<sse> {
  typedef int time_type;
  static constexpr bool is_path_integral = false;
  template<class Operator> static time_type time(const Operator&, int index) { return index; }
  static time_type end_time(std::size_t count) { return count; }
  static double flip_threshold(double, double, bool) { return 0; }
};

template<> struct worker_representation<path_integral> {
  typedef loop_config::time_t time_type;
  static constexpr bool is_path_integral = true;
  template<class Operator> static time_type time(const Operator& op, int) { return op.time(); }
  static time_type end_time(std::size_t) { return 1; }
  static double flip_threshold(double beta, double weight, bool field) {
    return field ? std::tanh(beta * weight) : 0;
  }
};

template<class Representation>
class loop_worker : public alps::parapack::mc_worker, private loop_config {
public:
  typedef Representation mc_type;
  typedef worker_representation<Representation> representation;
  typedef typename representation::time_type time_t;
  typedef looper::local_operator<mc_type, loop_graph_t, time_t> local_operator_t;
  typedef std::vector<local_operator_t> operator_string_t;
  typedef typename operator_string_t::iterator operator_iterator;

  typedef looper::union_find::node cluster_fragment_t;
  typedef looper::cluster_info cluster_info_t;

  typedef typename looper::estimator<measurement_set, mc_type, lattice_t, time_t>::type estimator_t;
  typedef double weight_parameter_type;

  loop_worker(alps::Parameters const& p);
  void init_observables(alps::Parameters const& params, alps::ObservableSet& obs);

  bool is_thermalized() const { return mcs.is_thermalized(); }
  double progress() const { return mcs.progress(); }

  void run(alps::ObservableSet& obs);

  // for exchange Monte Carlo
  void set_beta(double beta) { temperature.set_beta(beta); }
  double weight_parameter() const { return operators.size(); }
  static double log_weight(double gw, double beta) { return std::log(beta) * gw; }

  void save(alps::ODump& dp) const { dp << mcs << spins << operators; }
  void load(alps::IDump& dp) { dp >> mcs >> spins >> operators; }

protected:
  void build();
  void symmetrize();

  template<typename FIELD, typename SIGN, typename IMPROVE>
  void flip(alps::ObservableSet& obs);

private:
  // helpers
  lattice_t lattice;
  model_t model;

  // parameters
  looper::temperature temperature;
  double beta;
  bool use_improved_estimator;

  // configuration (checkpoint)
  looper::mc_steps mcs;
  std::vector<int> spins;
  std::vector<local_operator_t> operators;

  // observables
  double sign;
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

template<class Representation>
loop_worker<Representation>::loop_worker(alps::Parameters const& p)
  : alps::parapack::mc_worker(p), lattice(p), model(p, lattice, representation::is_path_integral),
    temperature(p), mcs(p) {

  if (temperature.annealing_steps() > mcs.thermalization())
    boost::throw_exception(std::invalid_argument("longer annealing steps than thermalization"));

  model.check_parameter(representation::is_path_integral && support_longitudinal_field, support_negative_sign);

  use_improved_estimator = (!model.has_field()) && (!p.defined("DISABLE_IMPROVED_ESTIMATOR"));
  if (!use_improved_estimator) std::cout << "WARNING: improved estimator is disabled\n";

  // configuration
  int nvs = num_sites(lattice.vg());
  spins.resize(nvs); std::fill(spins.begin(), spins.end(), 0 /* all up */);
  spins_c.resize(nvs);
  current.resize(nvs);
  perm.resize(max_virtual_sites(lattice));

  // initialize estimators
  estimator.initialize(p, lattice, model.is_signed(), use_improved_estimator);
}

template<class Representation>
void loop_worker<Representation>::init_observables(alps::Parameters const&, alps::ObservableSet& obs) {
  obs << make_observable(alps::SimpleRealObservable("Temperature"));
  obs << make_observable(alps::SimpleRealObservable("Inverse Temperature"));
  obs << make_observable(alps::SimpleRealObservable("Volume"));
  obs << make_observable(alps::SimpleRealObservable("Number of Sites"));
  obs << make_observable(alps::SimpleRealObservable("Number of Clusters"));
  if (model.is_signed()) {
    obs << alps::RealObservable("Sign");
    if (use_improved_estimator) {
      obs << alps::RealObservable("Weight of Zero-Meron Sector");
      obs << alps::RealObservable("Sign in Zero-Meron Sector");
    }
  }
  looper::energy_estimator::init_observables(obs, model.is_signed());
  estimator.init_observables(obs, model.is_signed());
}

template<class Representation>
void loop_worker<Representation>::run(alps::ObservableSet& obs) {
  // if (!mcs.can_work()) return;
  ++mcs;
  beta = 1.0 / temperature(mcs());

  build();
  symmetrize();

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

template<class Representation>
void loop_worker<Representation>::symmetrize() {
  int nvs = num_sites(lattice.vg());
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
        spins.begin() + offset, spins_c.begin() + offset, generator_01());
      for (int i = 0; i < s2; ++i) unify(fragments, offset+i, current[offset+perm[i]]);
    }
  }
}

//
// cluster flip
//

template<class Representation>
template<typename FIELD, typename SIGN, typename IMPROVE>
void loop_worker<Representation>::flip(alps::ObservableSet& obs) {
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
  int operator_index = 0;
  int negop = 0; // number of operators with negative weights
  BOOST_FOREACH(local_operator_t& op, operators) {
    time_t t = representation::time(op, operator_index++);
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
  }
  for (unsigned int s = 0; s < nvs; ++s) {
    weight.stop_top(current[s], representation::end_time(operators.size()), s, spins_c[s]);
    accum.stop_top(current[s], representation::end_time(operators.size()), s, spins_c[s]);
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
    to_flip[c] = ((2*uniform_01()-1) < representation::flip_threshold(beta, clusters[c].weight, FIELD()));
    if (SIGN() && IMPROVE() && (clusters[c].sign & 1) == 1) improved_sign = 0;
  }

  // improved measurement
  if (IMPROVE())
    estimator.improved_measurement(obs, lattice, beta, improved_sign, spins, operators,
      spins_c, fragments, coll);

  // flip operators & spins
  BOOST_FOREACH(local_operator_t& op, operators)
    if (to_flip[fragments[op.loop_0()].id()] ^ to_flip[fragments[op.loop_1()].id()]) op.flip();
  for (int s = 0; s < nvs; ++s) if (to_flip[fragments[s].id()]) spins[s] ^= 1;

  //
  // measurement
  //

  obs["Temperature"] << 1/beta;
  obs["Inverse Temperature"] << beta;
  obs["Volume"] << (double)lattice.volume();
  obs["Number of Sites"] << (double)num_sites(lattice.rg());
  obs["Number of Clusters"] << coll.num_clusters();

  // sign
  if (SIGN()) {
    if (IMPROVE()) {
      obs["Sign"] << improved_sign;
      if (alps::numeric::is_zero(improved_sign)) {
        obs["Weight of Zero-Meron Sector"] << 0.;
      } else {
        obs["Weight of Zero-Meron Sector"] << 1.;
        obs["Sign in Zero-Meron Sector"] << improved_sign;
      }
    } else {
      obs["Sign"] << sign;
    }
  }

  // energy
  double nop = coll.num_operators();
  double ene = model.energy_offset() - nop / beta;
  if (representation::is_path_integral && FIELD())
    for (unsigned int c = 0; c < clusters.size(); ++c)
      ene += (to_flip[c] ? -clusters[c].weight : clusters[c].weight);
  looper::energy_estimator::measurement(obs, lattice, beta, nop, sign, ene);

  // normal measurement
  estimator.normal_measurement(obs, lattice, beta, sign, spins, operators, spins_c);
}

} // namespace looper

#endif
