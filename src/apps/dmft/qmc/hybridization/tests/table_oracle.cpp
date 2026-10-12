// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
// Boundary regressions for retained CT-HYB tables. Affine Delta is reproduced
// exactly by linear interpolation; the quadratic K interpolation error is
// -(tau-left)*(right-tau). These are input-table and symmetry contracts, not
// Monte Carlo convergence or a complete retarded impurity reference.
#include <gtest/gtest.h>
#include <alps/hdf5/vector.hpp>
#include <alps/testing/temporary_directory.hpp>
#include "hybridization/hybfun.hpp"
#include "hybridization/hybretintfun.hpp"
#include <cmath>
#include <limits>
#include <sstream>

namespace {
alps::params parameters(const std::filesystem::path& filename, double beta, int cells) {
  alps::params p;
  p["BETA"]=beta; p["N_TAU"]=cells; p["N_ORBITALS"]=2;
  p["DELTA"]=filename.string(); p["DELTA_IN_HDF5"]=true;
  p["RET_INT_K"]=filename.string(); p["K_IN_HDF5"]=true;
  return p;
}
void write_tables(const std::filesystem::path& filename, double beta, int cells) {
  alps::hdf5::archive ar(filename.string(),"w");
  for (int f=0; f<2; ++f) {
    std::vector<double> delta;
    for (int i=0; i<=cells; ++i) delta.push_back(-(.2+f+.3*beta*i/cells));
    ar["/Delta_"+std::to_string(f)] << delta;
  }
  std::vector<double> kernel,derivative;
  for (int i=0; i<=cells; ++i) {
    const double tau=beta*i/cells;
    kernel.push_back(tau*(beta-tau)); derivative.push_back(beta-2.*tau);
  }
  ar["/Ret_int_K"] << kernel; ar["/Ret_int_Kp"] << derivative;
}
}

TEST(CthybTableOracle, AffineHybridizationHasDefinedEndpointsAndAntiperiodicInterior) {
  alps::testing::TemporaryDirectory tmp;
  for (double beta : {2.,5.}) for (int cells : {1,3,8,12}) {
    SCOPED_TRACE(beta);
    SCOPED_TRACE(cells);
    const auto file=tmp.path()/"tables.h5";
    write_tables(file,beta,cells);
    const hybfun delta(parameters(file,beta,cells));
    for (int f=0; f<2; ++f) {
      EXPECT_DOUBLE_EQ(delta.interpolate(0.,f),.2+f+.3*beta);
      EXPECT_DOUBLE_EQ(delta.interpolate(beta,f),.2+f);
      EXPECT_DOUBLE_EQ(delta.interpolate(-beta,f),-(.2+f+.3*beta));
      for (double fraction : {.01,.125,.3,.5,.9,.999}) {
        const double tau=beta*fraction;
        const double expected=.2+f+.3*(beta-tau);
        EXPECT_NEAR(delta.interpolate(tau,f),expected,1.e-13);
        EXPECT_NEAR(delta.interpolate(tau-beta,f),-expected,1.e-13);
      }
    }
    EXPECT_THROW(delta.interpolate(beta*1.01,0),std::invalid_argument);
    EXPECT_THROW(delta.interpolate(-beta*1.01,0),std::invalid_argument);
    EXPECT_THROW(delta.interpolate(std::numeric_limits<double>::quiet_NaN(),0),std::invalid_argument);
    EXPECT_THROW(delta.interpolate(std::numeric_limits<double>::infinity(),0),std::invalid_argument);
    EXPECT_THROW(delta.interpolate(0.,-1),std::out_of_range);
    EXPECT_THROW(delta.interpolate(0.,2),std::out_of_range);
  }
}

TEST(CthybTableOracle, RetardedKernelAndDerivativeHaveEvenAndOddEndpointContracts) {
  alps::testing::TemporaryDirectory tmp;
  for (double beta : {2.,5.}) for (int cells : {1,3,8,12}) {
    SCOPED_TRACE(beta);
    SCOPED_TRACE(cells);
    const auto file=tmp.path()/"tables.h5";
    write_tables(file,beta,cells);
    const ret_int_fun kernel(parameters(file,beta,cells));
    EXPECT_DOUBLE_EQ(kernel.interpolate(0.),0.);
    EXPECT_DOUBLE_EQ(kernel.interpolate(beta),0.);
    EXPECT_DOUBLE_EQ(kernel.interpolate(-beta),0.);
    EXPECT_DOUBLE_EQ(kernel.interpolate_deriv(0.),beta);
    EXPECT_DOUBLE_EQ(kernel.interpolate_deriv(beta),-beta);
    EXPECT_DOUBLE_EQ(kernel.interpolate_deriv(-beta),beta);
    for (int cell=0; cell<cells; ++cell) for (double fraction : {.25,.5,.75}) {
      const double left=beta*cell/cells,right=beta*(cell+1)/cells;
      const double tau=left+(right-left)*fraction;
      const double expected=tau*(beta-tau)-(tau-left)*(right-tau);
      EXPECT_NEAR(kernel.interpolate(tau),expected,1.e-12);
      EXPECT_NEAR(kernel.interpolate(-tau),expected,1.e-12);
      EXPECT_NEAR(kernel.interpolate_deriv(tau),beta-2.*tau,1.e-12);
      EXPECT_NEAR(kernel.interpolate_deriv(-tau),-(beta-2.*tau),1.e-12);
    }
    for (double invalid : {beta*1.01,-beta*1.01,
                           std::numeric_limits<double>::quiet_NaN(),
                           std::numeric_limits<double>::infinity()}) {
      EXPECT_THROW(kernel.interpolate(invalid),std::invalid_argument);
      EXPECT_THROW(kernel.interpolate_deriv(invalid),std::invalid_argument);
    }
  }
}

TEST(CthybTableOracle, ShortMeshDiagnosticsUseTheSuppliedStream) {
  alps::testing::TemporaryDirectory tmp;
  const auto file=tmp.path()/"tables.h5";
  write_tables(file,2.,1);
  const auto p=parameters(file,2.,1);
  const hybfun delta(p);
  const ret_int_fun kernel(p);
  std::ostringstream output;
  testing::internal::CaptureStdout();
  output << delta << kernel;
  const auto accidental_stdout=testing::internal::GetCapturedStdout();
  EXPECT_TRUE(accidental_stdout.empty());
  // Two mesh rows plus the repeated final row per table; never ten rows.
  EXPECT_NE(output.str().find("0 -0.2 -1.2"),std::string::npos);
  EXPECT_NE(output.str().find("1 -0.8 -1.8"),std::string::npos);
  EXPECT_NE(output.str().find("0 0 2"),std::string::npos);
  EXPECT_NE(output.str().find("1 0 -2"),std::string::npos);
  EXPECT_EQ(output.str().find("\n2 "),std::string::npos);
}
