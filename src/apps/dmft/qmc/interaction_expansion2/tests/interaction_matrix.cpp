// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
// Owned U-matrix input adapts pinned reconciliation90a36c0d4c64ee7308dcb0651264197ffdf36449
// qmc/U_matrix.h. Independent explicit density-interaction entries preserve the
// Gull/Werner/Fuchs/Troyer Hund convention, including the one-flavor legacy case.
#include <gtest/gtest.h>
#include <alps/testing/temporary_directory.hpp>
#include "U_matrix.h"

namespace {
alps::Parameters controls(int flavors,double interaction=4.,double hund=.5) {
  alps::Parameters p; p["FLAVORS"]=flavors; p["SITES"]=1; p["U"]=interaction; p["J"]=hund;
  return p;
}
}
TEST(CtintInteractionMatrix, GeneratedHundAndLegacySmallModelsHaveIndependentOwnedCopies) {
  const U_matrix matrix(controls(4));
  const std::vector<double> expected{0.,4.,2.5,3.,4.,0.,3.,2.5,2.5,3.,0.,4.,3.,2.5,4.,0.};
  EXPECT_EQ(matrix.data(),expected); EXPECT_EQ(matrix.n_nonzero(),12);
  EXPECT_DOUBLE_EQ(matrix.mu_shift(),4.75);
  auto copy=matrix; copy(0,1)=17.; EXPECT_DOUBLE_EQ(matrix(0,1),4.);
  const U_matrix single(controls(1));
  EXPECT_EQ(single.data(),(std::vector<double>{4.})); EXPECT_DOUBLE_EQ(single.mu_shift(),2.);
  const U_matrix two(controls(2));
  EXPECT_EQ(two.data(),(std::vector<double>{0.,4.,4.,0.})); EXPECT_DOUBLE_EQ(two.mu_shift(),2.);
  EXPECT_THROW((void)U_matrix(controls(3)),std::invalid_argument);
}
TEST(CtintInteractionMatrix, SparseFileAndCopiedMatrixSurviveSourceRemoval) {
  alps::testing::TemporaryDirectory tmp;
  const auto filename=tmp.path()/"interaction.txt";
  { std::ofstream file(filename); file << "0 1 2\n1 0 2\n0 2 -3\n2 0 -3\n"; }
  auto p=controls(5); p["U_MATRIX"]=filename.string();
  const U_matrix from_file(p);
  std::filesystem::remove(filename);
  const U_matrix captured(from_file);
  EXPECT_EQ(captured.data(),from_file.data()); EXPECT_EQ(captured.n_nonzero(),4);
  EXPECT_DOUBLE_EQ(captured.mu_shift(),-.5); EXPECT_DOUBLE_EQ(captured(4,0),0.);
  EXPECT_THROW((void)U_matrix(p),std::runtime_error);
}
TEST(CtintInteractionMatrix, MalformedIndicesSymmetryDiagonalAndExtentsFailExplicitly) {
  alps::testing::TemporaryDirectory tmp;
  const auto filename=tmp.path()/"interaction.txt";
  auto p=controls(4); p["U_MATRIX"]=filename.string();
  for (const char* text : {"-1 0 2\n","0 4 2\n","0 1 2\n1 0\n","0 1 2\ntrailing\n","0 1 2\n1 0 3\n","0 0 1\n"}) {
    SCOPED_TRACE(text); { std::ofstream file(filename); file << text; }
    EXPECT_THROW((void)U_matrix(p),std::invalid_argument);
  }
  EXPECT_THROW((void)U_matrix(controls(0)),std::invalid_argument);
  EXPECT_THROW((void)U_matrix(controls(-1)),std::invalid_argument);
  EXPECT_THROW((void)U_matrix(controls(std::numeric_limits<int>::max())),std::invalid_argument);
}
