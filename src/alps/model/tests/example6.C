/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2003-2004 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Axel Grzesik <axel@th.physik.uni-bonn.de>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include "basis_checks.hpp"

TEST(ModelBasis, DefaultSpinfulBosonStates) {
  const alps::site_basis<short> basis(model_test::descriptor("spinful boson"));
  model_test::expect_states(basis, {"N", "J", "Jz"},
                            {{0, 0, 0}, {1, 0, 0}, {1, 1, -1}, {1, 1, 0}, {1, 1, 1}});
}

TEST(ModelBasis, SpinfulBosonHonorsSpinAndOccupationCutoff) {
  alps::Parameters parameters;
  parameters["boson_spin"] = 2;
  parameters["NMax"] = 2;
  const alps::site_basis<short> basis(model_test::descriptor("spinful boson", parameters));
  ASSERT_EQ(basis.size(), 35u);
  model_test::expect_states(basis, {"N", "J", "Jz"}, model_test::spinful_states(2, 2));
}

TEST(ModelBasis, TJHasOnlyVacuumAndSingleParticleStates) {
  const alps::site_basis<short> basis(model_test::descriptor("t-J"));
  model_test::expect_states(basis, {"N", "S", "Sz"}, {{0, 0, 0}, {1, 0.5, -0.5}, {1, 0.5, 0.5}});
}

TEST(ModelBasis, AlternativeTJExcludesDoubleOccupancy) {
  const alps::site_basis<short> basis(model_test::descriptor("alternative t-J"));
  model_test::expect_states(basis, {"Nup", "Ndown"}, {{0, 0}, {0, 1}, {1, 0}});
}
