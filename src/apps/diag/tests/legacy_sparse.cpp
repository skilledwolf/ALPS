// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
// Legacy-file adaptation of reconciliation ef4b4ff36 tests/native_run.cpp.
// The four-site periodic Heisenberg spectrum is the independent oracle;
// compare eigenvalues and labels, not arbitrary degenerate eigenvectors.
#include "../sparsediag/sparsediag.h"
#include "legacy_reference.hpp"
TEST(EDLegacyRegression, SparseGroundStates) {
  alps::testing::TemporaryDirectory tmp;
  probe<SparseDiagMatrix<double>> sparse(input(tmp.path())); sparse.dostep();
  near(sparse.energies(),{-2,-1,-1,1,1});
}
TEST(EDLegacyRegression, ExhaustedLanczosThrowsInsteadOfPublishingSpectrum) {
  alps::testing::TemporaryDirectory tmp;
  probe<SparseDiagMatrix<double>> sparse(input(tmp.path(),{{"Sz_total","0"},{"MAX_ITERATIONS","1"}}));
  EXPECT_THROW(sparse.dostep(),std::runtime_error); EXPECT_TRUE(sparse.energies().empty());
}
