// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
// Sparse-pair proposal correction adapts reconciliation90a36c0d4c64ee7308dcb0651264197ffdf36449
// interaction_expansion2/model.cpp. These deterministic proposal partitions and
// exact first-vertex determinant ratios test retained-kernel detailed balance;
// they do not establish sampled interacting convergence or full multiband tails.
#include <gtest/gtest.h>
#include <alps/testing/temporary_directory.hpp>
#include "interaction_expansion2/interaction_expansion.hpp"
#include <fstream>
#include <map>
#include <limits>

namespace {
alps::params controls(int flavors=5) {
  alps::params p;
  p["BETA"]=2.; p["U"]=0.; p["J"]=0.; p["MU"]=0.; p["ALPHA"]=-.01;
  p["FLAVORS"]=flavors; p["SITES"]=1; p["NMATSUBARA"]=4; p["N"]=8;
  p["SWEEPS"]=64; p["THERMALIZATION"]=0; p["MEASUREMENT_PERIOD"]=1;
  p["RECALC_PERIOD"]=3; p["MAX_TIME"]=60.; p["SEED"]=42; p["ATOMIC"]=true;
  for (int f=0; f<flavors; ++f) { p["EPS_"+std::to_string(f)]=0.; p["EPSSQ_"+std::to_string(f)]=0.; }
  return p;
}
using pair=std::pair<spin_t,spin_t>;
class PairProbe : public MultiBandDensityHubbardInteractionExpansion {
public:
  using MultiBandDensityHubbardInteractionExpansion::MultiBandDensityHubbardInteractionExpansion;
  using MultiBandDensityHubbardInteractionExpansion::pair_at;
  const auto& pairs() const { return interaction_pairs_; }
  pair trial(double& insertion,double& removal) {
    insertion=try_add();
    if (vertices.empty()) throw std::logic_error("No trial vertex");
    const pair proposed{vertices.back().flavor1(),vertices.back().flavor2()};
    perform_add();
    sign=insertion<0. ? -1. : 1.;
    removal=try_remove(0);
    perform_remove(0);
    EXPECT_TRUE(vertices.empty());
    for (const auto& matrix : M) EXPECT_EQ(num_rows(matrix.matrix()),0);
    sign=1.;
    return proposed;
  }
};
std::filesystem::path matrix_file(const alps::testing::TemporaryDirectory& tmp,const std::map<pair,double>& interactions) {
  const auto path=tmp.path()/"interaction.txt";
  std::ofstream out(path);
  for (const auto& [edge,u] : interactions) out << edge.first << ' ' << edge.second << ' ' << u << '\n'
                                              << edge.second << ' ' << edge.first << ' ' << u << '\n';
  return path;
}
}

TEST(CtintPairOracle, EveryNonzeroUnorderedPairHasTheSameProposalInterval) {
  alps::testing::TemporaryDirectory tmp;
  // Unequal degrees (3,2,2,1) plus an isolated flavor. Row-conditioned retry
  // would over-propose edge0-3 and hang whenever the isolated flavor is drawn.
  const std::map<pair,double> interactions{{{0,1},2.},{{0,2},-3.},{{0,3},.7},{{1,2},1.2}};
  auto p=controls(); p["U_MATRIX"]=matrix_file(tmp,interactions).string();
  PairProbe worker(p,0);
  ASSERT_EQ(worker.pairs().size(),4);
  unsigned index=0;
  for (const auto& [edge,u] : interactions) {
    for (double fraction : {.01,.25,.5,.75,.99}) EXPECT_EQ(worker.pair_at((index+fraction)/4.),edge);
    ++index;
  }
  EXPECT_EQ(worker.pair_at(0.),interactions.begin()->first);
  EXPECT_THROW(worker.pair_at(1.),std::domain_error);
  EXPECT_THROW(worker.pair_at(std::numeric_limits<double>::quiet_NaN()),std::domain_error);
}

TEST(CtintPairOracle, ProductionInsertionRemovalRatiosAreReciprocalForSparseMixedSigns) {
  alps::testing::TemporaryDirectory tmp;
  const std::map<pair,double> interactions{{{0,1},2.},{{0,2},-3.},{{0,3},.7},{{1,2},1.2}};
  auto p=controls(); p["U_MATRIX"]=matrix_file(tmp,interactions).string();
  PairProbe worker(p,0);
  std::map<pair,unsigned> observed;
  for (unsigned trial=0; trial<128; ++trial) {
    double insertion,removal;
    const auto edge=worker.trial(insertion,removal);
    ASSERT_TRUE(interactions.count(edge)); ++observed[edge];
    // The half-filled zero-order determinant product is (-.51)*(.51).
    // Reciprocal edge proposal probability multiplies beta*U by4 pairs.
    const double expected=2.*interactions.at(edge)*4.*.2601;
    EXPECT_NEAR(insertion,expected,1.e-12);
    EXPECT_NEAR(insertion*removal,1.,1.e-12);
    EXPECT_EQ(insertion<0.,interactions.at(edge)<0.);
  }
  EXPECT_EQ(observed.size(),interactions.size());
}

TEST(CtintPairOracle, EmptyInteractionGraphAdvancesTheRetainedFreeWorkerWithoutVertices) {
  PairProbe worker(controls(4),0);
  EXPECT_TRUE(worker.pairs().empty());
  for (int step=0; step<64; ++step) { worker.update(); worker.measure(); }
  EXPECT_DOUBLE_EQ(worker.fraction_completed(),1.);
  const auto results=worker.collect_results();
  EXPECT_EQ(results["Sign"].count(),64); EXPECT_DOUBLE_EQ(results["Sign"].mean<double>(),1.);
  const auto density=results["densities"].mean<std::vector<double>>();
  ASSERT_EQ(density.size(),4);
  for (double n : density) EXPECT_DOUBLE_EQ(n,.5);
  for (double order : results["PertOrder"].mean<std::vector<double>>()) EXPECT_DOUBLE_EQ(order,0.);
}

TEST(CtintPairOracle, InvalidDensityInteractionsFailBeforeProposalSampling) {
  alps::testing::TemporaryDirectory tmp;
  auto p=controls();
  p["U_MATRIX"]=(tmp.path()/"interaction.txt").string();
  { std::ofstream out(tmp.path()/"interaction.txt"); out << "0 1 1\n1 0 2\n"; }
  EXPECT_THROW(PairProbe(p,0),std::invalid_argument);
  { std::ofstream out(tmp.path()/"interaction.txt"); out << "0 0 1\n"; }
  EXPECT_THROW(PairProbe(p,0),std::invalid_argument);
}
