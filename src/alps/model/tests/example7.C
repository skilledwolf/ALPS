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

TEST(ModelBasis, SpinfulBosonOccupationTwoStates) {
  alps::Parameters parameters;
  parameters["NMax"] = 2;
  const alps::site_basis<short> basis(model_test::descriptor("spinful boson", parameters));
  ASSERT_EQ(basis.size(), 14u);
  model_test::expect_states(basis, {"N", "J", "Jz"}, model_test::spinful_states(2, 1));
}

TEST(ModelBasis, AddingQuantumNumberDescriptorsExtendsOccupationRange) {
  alps::Parameters parameters;
  parameters["NMax"] = 2;
  auto descriptor = model_test::descriptor("spinful boson", parameters);
  const auto other = model_test::library().get_site_basis("spinful boson");
  ASSERT_EQ(descriptor.size(), other.size());
  for (std::size_t i = 0; i < descriptor.size(); ++i) descriptor[i] += other[i];
  const alps::site_basis<short> basis(descriptor);
  // Preserve the old example7 scenario: descriptor addition extends N to four.
  // This is enumeration of the combined ranges, not a tensor-product basis.
  ASSERT_EQ(basis.size(), 55u);
  model_test::expect_states(basis, {"N", "J", "Jz"}, model_test::spinful_states(4, 1));
}
