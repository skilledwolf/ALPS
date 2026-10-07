/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2003-2004 by Matthias Troyer <troyer@itp.phys.ethz.ch>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include "basis_checks.hpp"

TEST(ModelBasis, FermionEnumeratesBothOccupations) {
  const alps::site_basis<short> basis(model_test::descriptor("fermion"));
  model_test::expect_states(basis, {"Nup", "Ndown"}, {{0, 0}, {0, 1}, {1, 0}, {1, 1}});
}

TEST(ModelBasis, HardcoreBosonHasEmptyAndOccupiedStates) {
  const alps::site_basis<short> basis(model_test::descriptor("hardcore boson"));
  model_test::expect_states(basis, {"N"}, {{0}, {1}});
}

namespace {
void expect_spin_basis(double spin) {
  alps::Parameters parameters;
  parameters["local_spin"] = spin;
  const alps::site_basis<short> basis(model_test::descriptor("spin", parameters));
  model_test::States expected;
  // Exactly 2*S+1 projections, in ascending order, with fixed total spin S.
  for (int i = 0; i <= int(2 * spin); ++i) expected.push_back({spin, -spin + i});
  model_test::expect_states(basis, {"S", "Sz"}, expected);
}
} // namespace

TEST(ModelBasis, SpinHalfProjections) { expect_spin_basis(0.5); }
TEST(ModelBasis, SpinOneProjections) { expect_spin_basis(1); }
TEST(ModelBasis, SpinThreeHalvesProjections) { expect_spin_basis(1.5); }
TEST(ModelBasis, SpinTwoProjections) { expect_spin_basis(2); }
