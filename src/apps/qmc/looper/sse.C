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

#ifndef LOOPER_ONLY_PATH_INTEGRAL

namespace looper {

template<>
void loop_worker<sse>::build() {
  // initialize spin & operator information
  int nop = operators.size();
  std::copy(spins.begin(), spins.end(), spins_c.begin());
  std::swap(operators, operators_p); operators.resize(0);

  // initialize cluster information (setup cluster fragments)
  int nvs = num_sites(lattice.vg());
  fragments.resize(0); fragments.resize(nvs);
  for (int s = 0; s < nvs; ++s) current[s] = s;

  double bw = beta * model.graph_weight();
  bool try_gap = true;
  for (operator_iterator opi = operators_p.begin(); try_gap || opi != operators_p.end();) {

    // diagonal update & labeling
    if (try_gap) {
      if ((nop+1) * uniform_01() < bw) {
        loop_graph_t g = model.choose_graph(generator_01());
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
        if (bw * uniform_01() < nop) {
          --nop;
          ++opi;
          continue;
        } else {
          if (opi->is_site()) {
            opi->assign_graph(model.choose_diagonal(generator_01(), opi->loc(),
              spins_c[opi->pos()]));
          } else {
            opi->assign_graph(model.choose_diagonal(generator_01(), opi->loc(),
              spins_c[source(opi->pos(), lattice.vg())],
              spins_c[target(opi->pos(), lattice.vg())]));
          }
        }
      } else {
        if (opi->is_bond())
          opi->assign_graph(model.choose_offdiagonal(generator_01(), opi->loc(),
            spins_c[source(opi->pos(), lattice.vg())],
            spins_c[target(opi->pos(), lattice.vg())]));
      }
      operators.push_back(*opi);
      ++opi;
      try_gap = true;
    }

    operator_iterator oi = operators.end() - 1;
    if (oi->is_bond()) {
      int s0 = source(oi->pos(), lattice.vg());
      int s1 = target(oi->pos(), lattice.vg());
      if (oi->is_offdiagonal()) {
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
typedef looper::loop_worker<looper::sse> worker_type;

//
// dynamic registration to the factories
//

PARAPACK_REGISTER_ALGORITHM(worker_type, "loop; sse");
PARAPACK_REGISTER_ALGORITHM(alps::parapack::single_exchange_worker<worker_type>, "loop; sse; exchange");

} // end namespace

#endif // LOOPER_ONLY_PATH_INTEGRAL
