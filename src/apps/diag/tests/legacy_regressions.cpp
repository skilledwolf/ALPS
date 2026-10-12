// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
// Legacy-file adaptation of reconciliation ef4b4ff36 tests/native_run.cpp.
// The four-site periodic Heisenberg spectrum is the independent oracle;
// compare eigenvalues and labels, not arbitrary degenerate eigenvectors.
#include "../fulldiag/fulldiag.h"
#include "legacy_reference.hpp"
TEST(EDLegacyRegression, CompleteRealSpectrum) {
  alps::testing::TemporaryDirectory tmp; const auto path=input(tmp.path());
  probe<FullDiagMatrix<double>> full(path); full.dostep();
  near(full.energies(),{-2,-1,-1,-1,0,0,0,0,0,0,0,1,1,1,1,1});
}
TEST(EDLegacyRegression, FixedQuantumNumberAndMomentumKeepTheirLabels) {
  alps::testing::TemporaryDirectory tmp;
  probe<FullDiagMatrix<std::complex<double>>> full(input(tmp.path(),
    {{"Sz_total","0"},{"TRANSLATION_SYMMETRY","1"},{"TOTAL_MOMENTUM","Pi"}}));
  full.dostep(); near(full.energies(),{-1,0}); ASSERT_EQ(full.labels().size(),1u);
  const auto& labels=full.labels()[0];
  EXPECT_NE(std::find(labels.begin(),labels.end(),std::make_pair(std::string("Sz"),std::string("0"))),labels.end());
  EXPECT_NE(std::find(labels.begin(),labels.end(),std::make_pair(std::string("TOTAL_MOMENTUM"),std::string("Pi"))),labels.end());
}
