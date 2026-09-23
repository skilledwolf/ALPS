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

#include "loop_worker.h"

namespace looper {

template<>
void loop_worker<path_integral>::build() {
  // initialize spin & operator information
  std::copy(spins.begin(), spins.end(), spins_c.begin());
  std::swap(operators, operators_p); operators.resize(0);

  // initialize cluster information (setup cluster fragments)
  int nvs = num_sites(lattice.vg());
  fragments.resize(0); fragments.resize(nvs);
  for (int s = 0; s < nvs; ++s) current[s] = s;

  boost::variate_generator<engine_type&, boost::exponential_distribution<> >
    r_time(engine(), boost::exponential_distribution<>(beta * model.graph_weight()));
  double t = r_time();
  for (operator_iterator opi = operators_p.begin(); t < 1 || opi != operators_p.end();) {

    // diagonal update & labeling
    if (opi == operators_p.end() || t < opi->time()) {
      loop_graph_t g = model.choose_graph(generator_01());
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

    operator_iterator oi = operators.end() - 1;
    if (oi->is_bond()) {
      int s0 = source(oi->pos(), lattice.vg());
      int s1 = target(oi->pos(), lattice.vg());
      if (oi->is_offdiagonal()) {
        oi->assign_graph(model.choose_offdiagonal(generator_01(), oi->loc(),
          spins_c[s0], spins_c[s1]));
        spins_c[s0] ^= 1;
        spins_c[s1] ^= 1;
      }
      boost::tie(current[s0], current[s1], oi->loop0, oi->loop1) =
        reconnect(fragments, oi->graph(), current[s0], current[s1]);
    } else {
      int s = oi->pos();
      if (oi->is_offdiagonal()) spins_c[s] ^= 1;
      boost::tie(current[s], oi->loop0, oi->loop1) = reconnect(fragments, oi->graph(), current[s]);
    }
  }

}

} // namespace looper

namespace {
typedef looper::loop_worker<looper::path_integral> worker_type;

typedef looper::evaluator<loop_config::measurement_set> loop_evaluator;

//
// dynamic registration to the factories
//

PARAPACK_REGISTER_ALGORITHM(worker_type, "loop");
PARAPACK_REGISTER_ALGORITHM(worker_type, "loop; path integral");
PARAPACK_REGISTER_ALGORITHM(alps::parapack::single_exchange_worker<worker_type>, "loop; exchange");
PARAPACK_REGISTER_ALGORITHM(alps::parapack::single_exchange_worker<worker_type>, "loop; path integral; exchange");
PARAPACK_REGISTER_EVALUATOR(loop_evaluator, "loop");

} // end namespace
