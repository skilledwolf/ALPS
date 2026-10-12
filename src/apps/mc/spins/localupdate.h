/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1999-2006 by Matthias Troyer <troyer@comp-phys.org>,
*                            Fabian Stoeckli <fabstoec@student.ethz.ch>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#ifndef ALPS_APPLICATIONS_MC_SPIN_LOCALUPDATE_H_
#define ALPS_APPLICATIONS_MC_SPIN_LOCALUPDATE_H_

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/random/uniform_int.hpp>
#include <vector>
#include <set>

#include "potts.h"
#include "ising.h"
#include "on.h"
#include "matrices.h"

// Full matrices are oriented by the original XML bond, as in the independent
// Hamiltonian used by reconciliation90a36c0d4c64ee7308dcb0651264197ffdf36449,
// src/apps/mc/spins/spinsim.h::edge_energy/local_update. Boost's undirected
// out_edges instead puts the currently updated vertex first. At the original
// target, x^T J y therefore requires J^T in the retained first-spin proposal.
// This is a bespoke correction of the retained Troyer/Stoeckli local kernel;
// scalar/diagonal matrices and their trajectories are unchanged.
template<class MAT> struct spinmc_oriented_matrix { static constexpr bool value = false; };
template<class T, int N> struct spinmc_oriented_matrix<SquareMatrix<T,N>> { static constexpr bool value = true; };

template<class MAT, class Graph>
auto spinmc_local_bond_sources(const Graph &graph) {
  using vertex = typename boost::graph_traits<Graph>::vertex_descriptor;
  std::vector<vertex> sources;
  if constexpr (spinmc_oriented_matrix<MAT>::value) {
    sources.resize(boost::num_edges(graph));
    auto edges = boost::edges(graph);
    for (auto edge = edges.first; edge != edges.second; ++edge)
      sources.at(boost::get(boost::edge_index,graph,*edge)) = boost::source(*edge,graph);
  }
  return sources;
}

template<class MAT, class Graph, class Edge, class Sources>
MAT spinmc_local_bond_matrix(MAT matrix, const Graph &graph, const Edge &edge, const Sources &sources) {
  if constexpr (spinmc_oriented_matrix<MAT>::value) {
    if (boost::source(edge,graph) != sources.at(boost::get(boost::edge_index,graph,edge))) {
      const auto original = matrix;
      for (int i = 0; i < MAT::dim; ++i)
        for (int j = 0; j < MAT::dim; ++j)
          matrix.setElement(i,j,original.getElement(j,i));
    }
  }
  return matrix;
}

// A self-bond contributes -s^T J s once. Its antisymmetric part vanishes;
// reuse the retained onsite reflection difference for its symmetric part.
template<class MAT> MAT spinmc_self_bond_matrix(MAT matrix) {
  if constexpr (spinmc_oriented_matrix<MAT>::value) {
    const auto original=matrix;
    for(int i=0;i<MAT::dim;++i) for(int j=0;j<MAT::dim;++j)
      matrix.setElement(i,j,.5*original.getElement(i,j)+.5*original.getElement(j,i));
  }
  return matrix;
}

template <class Graph, class MomentMap, class RNG, class CouplingMap, class SpinfactorMap>
void local_update(const Graph& graph, MomentMap& moment,
                  double beta, RNG& rng,  const CouplingMap& coupling,
                  const SpinfactorMap& spinfactor, double g,
                  TinyVector<double,CouplingMap::value_type::dim>
                    h = TinyVector<double,CouplingMap::value_type::dim>(0.0))

{
  // some typedefs for clarity
  typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename boost::property_traits<MomentMap>::value_type moment_type;
  typedef typename moment_type::update_type update_type;
  typedef typename CouplingMap::value_type MAT;
  const auto bond_sources = spinmc_local_bond_sources<MAT>(graph);

  // a local update consists of on averaqe one update per site
  for (int steps=0; steps<boost::num_vertices(graph);++steps) {
    // step i)
    vertex_descriptor vertex=*(boost::vertices(graph).first 
                              + int (rng()*boost::num_vertices(graph)));

    // step ii)
    update_type update=moment[vertex].random_update(rng);

    // step iii)
    double delta_E=0;
    std::set<typename boost::graph_traits<Graph>::edge_descriptor> self_bonds;
    typename boost::graph_traits<Graph>::out_edge_iterator edge,end;
    for (boost::tie(edge,end)=boost::out_edges(vertex,graph);edge!=end;++edge) {
      double fact = -spinfactor[boost::source(*edge,graph)]*spinfactor[boost::target(*edge,graph)];
      MAT J = spinmc_local_bond_matrix(MAT(coupling[*edge]),graph,*edge,bond_sources);
      J = J * fact;
      if (boost::source(*edge,graph)==boost::target(*edge,graph)) {
        if(self_bonds.insert(*edge).second) {
          J=spinmc_self_bond_matrix(J);
          delta_E+=onsite_energy_change(moment[vertex],J,update);
        }
        continue;
      }
          delta_E += bond_energy_change(moment[boost::source(*edge,graph)],
                               moment[boost::target(*edge,graph)],J,update);
    }
    delta_E += spinfactor[vertex]*site_energy_change(moment[vertex], h, update)  * g;
    
    // step iv)
    if (delta_E<0 || rng() < exp(-beta*delta_E)) 
      moment[vertex].update(update);
  }
}

// generalized version, with onsite_energy_terms
template <class Graph, class MomentMap, class RNG, class CouplingMap, 
            class SelfinteractionMap, class SpinfactorMap>
void local_update_self(const Graph& graph, MomentMap& moment,
                  double beta, RNG& rng,  const CouplingMap& coupling, 
                  const SelfinteractionMap& selfinteraction, 
                  const SpinfactorMap& spinfactor, double g,
                  TinyVector<double,CouplingMap::value_type::dim> 
                    h = TinyVector<double,CouplingMap::value_type::dim>(0.0))
{
  // some typedefs for clarity
  typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename boost::property_traits<MomentMap>::value_type moment_type;
  typedef typename moment_type::update_type update_type;
  typedef typename CouplingMap::value_type MAT;
  const auto bond_sources = spinmc_local_bond_sources<MAT>(graph);

  // a local update consists of on averaqe one update per site
  for (int steps=0; steps<boost::num_vertices(graph);++steps) {
    // step i)
    vertex_descriptor vertex=*(boost::vertices(graph).first + 
                                 int(rng()*boost::num_vertices(graph)));
    
    // step ii)
    update_type update = moment[vertex].random_update(rng);
  
    // step iii)
    double delta_E=0;
    std::set<typename boost::graph_traits<Graph>::edge_descriptor> self_bonds;
    typename boost::graph_traits<Graph>::out_edge_iterator edge,end;
    for (boost::tie(edge,end)=boost::out_edges(vertex,graph);edge!=end;++edge) {
      double fact = -spinfactor[boost::source(*edge,graph)]*spinfactor[boost::target(*edge,graph)];
      MAT J = spinmc_local_bond_matrix(MAT(coupling[*edge]),graph,*edge,bond_sources);
      J = J * fact;
      if (boost::source(*edge,graph)==boost::target(*edge,graph)) {
        if(self_bonds.insert(*edge).second) {
          J=spinmc_self_bond_matrix(J);
          delta_E+=onsite_energy_change(moment[vertex],J,update);
        }
        continue;
      }
      delta_E += bond_energy_change(moment[boost::source(*edge,graph)],
                                moment[boost::target(*edge,graph)],J,update);
    }
    delta_E += spinfactor[vertex]*site_energy_change(moment[vertex], h, update) * g;
    
    delta_E += onsite_energy_change(moment[vertex], selfinteraction[vertex],update);
    
    // step iv)
    if (delta_E<0 || rng() < exp(-beta*delta_E))
      moment[vertex].update(update); 
  }
}

#endif
