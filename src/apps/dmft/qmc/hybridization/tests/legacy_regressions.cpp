// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
// Legacy-interface regressions for fixes extracted from reconciliation
// 67fb28558 (warmup), b9de725c7 (G2 buffers), and 2346d1b76 (absent table).
#include <gtest/gtest.h>
#include <alps/testing/temporary_directory.hpp>
#include "hybridization/hyb.hpp"
#include <fstream>

TEST(CthybLegacyRegression, EveryWarmupBlockIsExcluded) {
  alps::testing::TemporaryDirectory tmp;
  const auto delta=tmp.path()/"delta.dat";
  { std::ofstream out(delta); for(int i=0;i<=4;++i) out<<i<<" -.2 -.2\n"; }
  for (int warmup : {0,2}) {
    alps::params p;
    p["BETA"]=2.; p["U"]=0.; p["MU"]=0.; p["N_ORBITALS"]=2;
    p["N_TAU"]=4; p["N_MEAS"]=3; p["SWEEPS"]=4; p["THERMALIZATION"]=warmup;
    p["SEED"]=42; p["MAX_TIME"]=60; p["DELTA"]=delta.string(); p["MEASURE_time"]=false;
    hybridization worker(p,0);
    for(int sweep=1;sweep<=warmup+4;++sweep) {
      worker.update(); worker.measure();
      EXPECT_EQ(alps::mcresult(worker.get_measurements()["Sign"]).count(),
                static_cast<unsigned>(std::max(0,sweep-warmup)));
    }
  }
}
TEST(CthybLegacyRegression, G2OnlyAllowsAnEmptyF2Buffer) {
  alps::params p; p["BETA"]=2.; p["MEASURE_g2w"]=true;
  hybmatrix matrix(p);
  std::vector<std::complex<double>> g2(4,{9.,1.}), f2;
  matrix.measure_G2w(g2,f2,2,2,{});
  EXPECT_TRUE(f2.empty());
  for(const auto value:g2) EXPECT_EQ(value,std::complex<double>());
}
TEST(CthybLegacyRegression, AbsentRetardedTableHasZeroValues) {
  alps::params p; p["BETA"]=2.; p["N_TAU"]=4;
  ret_int_fun table(p);
  for(double t : {0.,.5,2.}) {
    EXPECT_DOUBLE_EQ(table.interpolate(t),0.);
    EXPECT_DOUBLE_EQ(table.interpolate_deriv(t),0.);
  }
}
