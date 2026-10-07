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

TEST(ModelBasis, SingleQuantumNumberHardcoreBoson) {
  const alps::site_basis<short, alps::single_qn_site_state<short>> basis(
      model_test::descriptor("hardcore boson"));
  ASSERT_NO_FATAL_FAILURE(model_test::expect_states(basis, {"N"}, {{0}, {1}}));
  EXPECT_FALSE(alps::is_fermionic(basis, 0));
  EXPECT_FALSE(alps::is_fermionic(basis, 1));
}

TEST(ModelBasis, SingleQuantumNumberSpinlessFermion) {
  const alps::site_basis<short, alps::single_qn_site_state<short>> basis(
      model_test::descriptor("spinless fermion"));
  ASSERT_NO_FATAL_FAILURE(model_test::expect_states(basis, {"N"}, {{0}, {1}}));
  EXPECT_FALSE(alps::is_fermionic(basis, 0));
  EXPECT_TRUE(alps::is_fermionic(basis, 1));
}
