// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
// Regression oracles for the retained Gull/Werner/Fuchs/Troyer Hirsch-Fye
// implementation. Corrections were identified in reconciliation 90a36c0d4;
// the small linear systems and proposal probabilities below are analytic.
#include "../hirschfyeaux.h"
#include <gtest/gtest.h>
#include <array>
#include <limits>

namespace {
struct proposal {
  std::array<double, 2> draws;
  unsigned next = 0;
  double operator()() { return draws.at(next++); }
};
}

TEST(HirschFyeNumerics, NegativeDeterminantUsesAbsoluteWeightHeatBath) {
  // exp(2 lambda)=2, r_up=1, r_down=-1/2. The signed target is sampled
  // by its absolute weight, so P(accept)=1/3, not abs((-1/2)/(1-1/2))=1.
  dense_matrix up(1, 1), down(1, 1);
  up(0, 0)=1.; down(0, 0)=2.5;
  std::vector<int> spins{1};
  int sign=1;
  proposal reject{{0., .5}};
  update_single_spin(reject, up, down, spins, std::log(2.)/2., sign);
  EXPECT_EQ(spins[0], 1);
  EXPECT_EQ(sign, 1);
  EXPECT_DOUBLE_EQ(down(0, 0), 2.5);
  proposal accept{{0., .25}};
  update_single_spin(accept, up, down, spins, std::log(2.)/2., sign);
  EXPECT_EQ(spins[0], -1);
  EXPECT_EQ(sign, -1);
  EXPECT_NEAR(up(0, 0), 1., 1e-14);
  EXPECT_NEAR(down(0, 0), -5., 1e-13);
}

TEST(HirschFyeNumerics, ReconstructionSolvesPivotedSignedSystem) {
  // With exp(lambda)=2 and spins=(1,1), A=2I-G0=[[0,-1],[-3,1]].
  // det(A)=-3; A^-1 G0=[[-5/3,-2/3],[-2,-1]]. The first pivot swaps rows.
  dense_matrix bare(2, 2), result;
  bare(0,0)=2.; bare(0,1)=1.; bare(1,0)=3.; bare(1,1)=1.;
  EXPECT_EQ(update_from_zero(result, bare, {1,1}, std::log(2.)), -1);
  ASSERT_EQ(result.size1(), 2u);
  EXPECT_NEAR(result(0,0), -5./3., 1e-13);
  EXPECT_NEAR(result(0,1), -2./3., 1e-13);
  EXPECT_NEAR(result(1,0), -2., 1e-13);
  EXPECT_NEAR(result(1,1), -1., 1e-13);
}

TEST(HirschFyeNumerics, FailedReconstructionPreservesDestination) {
  dense_matrix bare(1,1), result(1,1);
  bare(0,0)=2.; result(0,0)=17.; // A=2-G0=0, singular.
  EXPECT_THROW(update_from_zero(result, bare, {1}, std::log(2.)), std::runtime_error);
  EXPECT_DOUBLE_EQ(result(0,0), 17.);
  bare(0,0)=std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(update_from_zero(result, bare, {1}, .1), std::invalid_argument);
  EXPECT_DOUBLE_EQ(result(0,0), 17.);
}

TEST(HirschFyeNumerics, RankOneFlipAgreesWithFullReconstruction) {
  dense_matrix bare(2,2), up, down, reference_up, reference_down;
  bare(0,0)=.5; bare(0,1)=.2; bare(1,0)=-.2; bare(1,1)=.5;
  std::vector<int> spins{1,-1};
  const double lambda=.4;
  int sign=update_from_zero(up, bare, spins, lambda)
          *update_from_zero(down, bare, spins, -lambda);
  proposal accept{{0., 0.}};
  update_single_spin(accept, up, down, spins, lambda, sign);
  ASSERT_EQ(spins[0], -1);
  const int reference_sign=update_from_zero(reference_up, bare, spins, lambda)
                          *update_from_zero(reference_down, bare, spins, -lambda);
  EXPECT_EQ(sign, reference_sign);
  for(unsigned i=0;i<2;++i) for(unsigned j=0;j<2;++j) {
    EXPECT_NEAR(up(i,j), reference_up(i,j), 2e-14);
    EXPECT_NEAR(down(i,j), reference_down(i,j), 2e-14);
  }
}

TEST(HirschFyeNumerics, NormCoversOtherColumnsAndNonfiniteEntries) {
  dense_matrix matrix(3,3,0.);
  matrix(2,2)=-7.;
  EXPECT_DOUBLE_EQ(norm_max(matrix), 7.);
  matrix(0,1)=std::numeric_limits<double>::quiet_NaN();
  EXPECT_EQ(norm_max(matrix), std::numeric_limits<double>::infinity());
}

TEST(HirschFyeNumerics, HubbardStratonovichSmallCouplingAndExponentBounds) {
  EXPECT_DOUBLE_EQ(hirschfye_lambda(1.,0.,1),0.);
  EXPECT_NEAR(hirschfye_lambda(1.,1e-20,1),1e-10,1e-20);
  EXPECT_NEAR(hirschfye_lambda(4.,2.,8),std::acosh(std::exp(.5)),1e-14);
  EXPECT_THROW(hirschfye_lambda(1.,-1.,8),std::invalid_argument);
  EXPECT_THROW(hirschfye_lambda(1.,1.,0),std::invalid_argument);
  EXPECT_THROW(hirschfye_lambda(1000.,1.,1),std::invalid_argument);
}
