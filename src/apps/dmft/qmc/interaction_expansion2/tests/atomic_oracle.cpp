// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
// The free-atom oracle is adapted from ALPS reconciliation commit
// 90a36c0d4c64ee7308dcb0651264197ffdf36449,
// src/apps/dmft/qmc/interaction_expansion2/tests/run_config.cpp and
// tests/multiband.cpp. This fixture exercises the retained two-flavor NGS
// solver and its production postprocessor, without migrating either API.
// At U=MU=0 the four occupation states are equiprobable, so n_f=1/2,
// n_f^2=1/2, n_up*n_down=1/4 and G(iw)=1/iw, G(tau)=-1/2.

#include <gtest/gtest.h>
#include <alps/hdf5/complex.hpp>
#include <alps/hdf5/vector.hpp>
#include <alps/testing/temporary_directory.hpp>
#include "interaction_expansion2/interaction_expansion.hpp"
#include <cmath>
#include <complex>
#include <filesystem>
#include <vector>

void compute_greens_functions(const alps::results_type<HubbardInteractionExpansion>::type&,
                             const alps::params&, const std::string&);

namespace {
// The retained postprocessor writes matrix_size in cwd. Confine that diagnostic
// alongside the explicitly named input/output archives, restoring cwd on error.
class ScopedDirectory {
public:
  explicit ScopedDirectory(const std::filesystem::path& directory)
    : previous_(std::filesystem::current_path()) {
    std::filesystem::current_path(directory);
  }
  ~ScopedDirectory() {
    std::error_code ignored;
    std::filesystem::current_path(previous_, ignored);
  }
private:
  std::filesystem::path previous_;
};

alps::params parameters(const std::string& input, double beta, bool atomic) {
  alps::params p;
  p["BETA"]=beta; p["U"]=0.; p["MU"]=0.; p["ALPHA"]=-.01;
  p["FLAVORS"]=2; p["SITES"]=1; p["NMATSUBARA"]=8; p["N"]=8;
  p["SWEEPS"]=64; p["THERMALIZATION"]=0; p["MEASUREMENT_PERIOD"]=1;
  p["RECALC_PERIOD"]=3; p["MAX_TIME"]=60; p["SEED"]=42;
  p["EPS_0"]=0.; p["EPS_1"]=0.; p["EPSSQ_0"]=0.; p["EPSSQ_1"]=0.;
  p["INFILE"]=input;
  if (atomic) p["ATOMIC"]=true;
  return p;
}
}

TEST(CtintAtomicOracle, OccupationsAndGreenFunctionsMatchTheFreeAtom) {
  alps::testing::TemporaryDirectory directory;
  ScopedDirectory cwd(directory.path());
  for (bool atomic : {false,true}) {
    for (double beta : {2.,5.}) {
      SCOPED_TRACE(atomic ? "atomic initialization" : "HDF5 frequency initialization");
      SCOPED_TRACE(beta);
      std::vector<std::complex<double>> bare(8);
      for (std::size_t w=0; w<bare.size(); ++w)
        bare[w]={0.,-beta/((2.*w+1.)*std::acos(-1.))};
      const auto input=(directory.path()/"g0.h5").string();
      {
        alps::hdf5::archive archive(input,"w");
        archive["/G0_0"] << bare;
        archive["/G0_1"] << bare;
      }
      const auto p=parameters(input,beta,atomic);
      HubbardInteractionExpansion simulation(p,0);
      // Explicit finite steps isolate the numerical estimator from wall-clock,
      // MPI and legacy stopping policy. Every insertion has zero weight at U=0.
      for (int step=0; step<64; ++step) {
        simulation.update();
        simulation.measure();
      }
      const auto results=simulation.collect_results();
      EXPECT_EQ(results["Sign"].count(),64);
      EXPECT_DOUBLE_EQ(results["Sign"].mean<double>(),1.);
      EXPECT_NEAR(results["density_correlation"].mean<double>(),.25,1.e-12);
      const auto density=results["densities"].mean<std::vector<double>>();
      ASSERT_EQ(density.size(),2);
      for (double value : density) EXPECT_NEAR(value,.5,1.e-12);
      const auto pairs=results["n_i n_j"].mean<std::vector<double>>();
      ASSERT_EQ(pairs.size(),4);
      for (std::size_t i=0; i<pairs.size(); ++i)
        EXPECT_NEAR(pairs[i],i==0 || i==3 ? .5 : .25,1.e-12);
      const auto order=results["PertOrder"].mean<std::vector<double>>();
      ASSERT_EQ(order.size(),2);
      for (double value : order) EXPECT_DOUBLE_EQ(value,0.);
      for (int flavor=0; flavor<2; ++flavor) {
        const auto suffix=std::to_string(flavor)+"_0_0";
        for (const auto* prefix : {"Wk_real_","Wk_imag_"}) {
          const auto values=results[std::string(prefix)+suffix].mean<std::vector<double>>();
          ASSERT_EQ(values.size(),bare.size());
          for (double value : values) EXPECT_DOUBLE_EQ(value,0.);
        }
      }
      const auto output=(directory.path()/"result.h5").string();
      compute_greens_functions(results,p,output);
      alps::hdf5::archive archive(output,"r");
      for (int flavor=0; flavor<2; ++flavor) {
        const auto suffix=std::to_string(flavor)+"/mean/value";
        std::vector<std::complex<double>> frequency;
        std::vector<double> time;
        archive["/G_omega/"+suffix] >> frequency;
        archive["/G_tau/"+suffix] >> time;
        ASSERT_EQ(frequency.size(),8);
        ASSERT_EQ(time.size(),9);
        for (std::size_t w=0; w<frequency.size(); ++w)
          EXPECT_NEAR(std::abs(frequency[w]-bare[w]),0.,1.e-12);
        for (double value : time) EXPECT_NEAR(value,-.5,1.e-12);
      }
    }
  }
}
