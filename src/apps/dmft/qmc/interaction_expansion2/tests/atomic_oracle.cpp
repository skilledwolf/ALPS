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

class FrequencyProbe : public HubbardInteractionExpansion {
public:
  using HubbardInteractionExpansion::HubbardInteractionExpansion;
  void measure_one_vertex() {
    // At U=4 and beta=2 the first-vertex proposal weight is +2.0808,
    // independent of time or symmetrized alpha: -beta*U*(-.51)*(.51).
    const double proposal=try_add();
    EXPECT_NEAR(proposal,2.0808,1.e-12);
    perform_add();
    sign=proposal<0 ? -1. : 1.;
    measure_observables();
  }
};

class RebuildProbe : public HubbardInteractionExpansion {
public:
  using HubbardInteractionExpansion::HubbardInteractionExpansion;
  void check_rebuild(double interaction) {
    const double proposal=try_add();
    // Exact first-vertex determinant product: (-.51)*(.51), independent of
    // time and symmetrized alpha. Both signs are unconditionally accepted.
    ASSERT_NEAR(proposal,.5202*interaction,1.e-12);
    ASSERT_GT(std::abs(proposal),1.);
    perform_add();
    sign=proposal<0 ? -1. : 1.;
    const double accepted_sign=sign;
    const auto before=M;
    ASSERT_EQ(vertices.size(),1);
    measure_observables();
    const auto first=collect_results();
    reset_perturbation_series();
    ASSERT_EQ(vertices.size(),1);
    EXPECT_EQ(sign,accepted_sign);
    for (unsigned flavor=0; flavor<2; ++flavor) {
      ASSERT_EQ(num_rows(M[flavor].matrix()),1);
      EXPECT_NEAR(std::abs(M[flavor].matrix()(0,0)),1./.51,1.e-12);
      EXPECT_NEAR(M[flavor].matrix()(0,0),before[flavor].matrix()(0,0),1.e-12);
    }
    measure_observables();
    const auto after=collect_results();
    EXPECT_EQ(after["Sign"].count(),2);
    EXPECT_EQ(after["Sign"].mean<double>(),accepted_sign);
    const auto first_density=first["densities"].mean<std::vector<double>>();
    const auto after_density=after["densities"].mean<std::vector<double>>();
    ASSERT_EQ(first_density.size(),2);
    ASSERT_EQ(after_density.size(),2);
    for (unsigned flavor=0; flavor<2; ++flavor)
      EXPECT_NEAR(after_density[flavor],first_density[flavor],1.e-12);
  }
};
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

TEST(CtintAtomicOracle, ExcludesEveryWarmupConfigurationAndReportsFiniteCompletion) {
  // Adapted boundary from pinned reconciliation tests/multiband.cpp. This
  // retained two-flavor oracle counts elementary updates, without MPI/time.
  auto p=parameters("unused",2.,true);
  p["THERMALIZATION"]=3;
  p["SWEEPS"]=4;
  HubbardInteractionExpansion simulation(p,0);
  EXPECT_FALSE(simulation.is_thermalized());
  EXPECT_DOUBLE_EQ(simulation.fraction_completed(),0.);
  simulation.measure();
  EXPECT_EQ(simulation.collect_results()["Sign"].count(),0);
  for (int warmup=0; warmup<3; ++warmup) {
    simulation.update(); simulation.measure();
    EXPECT_FALSE(simulation.is_thermalized());
    EXPECT_DOUBLE_EQ(simulation.fraction_completed(),0.);
    EXPECT_EQ(simulation.collect_results()["Sign"].count(),0);
  }
  for (int physical=1; physical<=4; ++physical) {
    simulation.update(); simulation.measure();
    EXPECT_TRUE(simulation.is_thermalized());
    EXPECT_DOUBLE_EQ(simulation.fraction_completed(),physical/4.);
    const auto results=simulation.collect_results();
    EXPECT_EQ(results["Sign"].count(),physical);
    EXPECT_EQ(results["densities"].count(),physical);
    EXPECT_DOUBLE_EQ(results["Sign"].mean<double>(),1.);
  }
}

TEST(CtintAtomicOracle, OrdinaryRunStopsAtTheRequestedPostWarmupWork) {
  auto p=parameters("unused",2.,true);
  p["THERMALIZATION"]=3;
  p["SWEEPS"]=4;
  HubbardInteractionExpansion simulation(p,0);
  unsigned polls=0;
  // Bound a regression of the old never-completing path without wall-clock
  // waiting. Correct completion occurs before this independent safety limit.
  EXPECT_TRUE(simulation.run([&] { return ++polls>16; }));
  EXPECT_LT(polls,16u);
  EXPECT_DOUBLE_EQ(simulation.fraction_completed(),1.);
  EXPECT_EQ(simulation.collect_results()["Sign"].count(),4);
}

TEST(CtintAtomicOracle, MeasurementBlocksThatCrossWarmupCountTheFinalConfiguration) {
  auto p=parameters("unused",2.,true);
  p["THERMALIZATION"]=4;
  p["SWEEPS"]=5;
  p["MEASUREMENT_PERIOD"]=3;
  HubbardInteractionExpansion simulation(p,0);
  simulation.update(); simulation.measure(); // step=3: warmup
  EXPECT_EQ(simulation.collect_results()["Sign"].count(),0);
  EXPECT_DOUBLE_EQ(simulation.fraction_completed(),0.);
  simulation.update(); simulation.measure(); // step=6: one final configuration
  EXPECT_EQ(simulation.collect_results()["Sign"].count(),1);
  EXPECT_DOUBLE_EQ(simulation.fraction_completed(),.4);
  simulation.update(); simulation.measure(); // step=9: requested work completed
  EXPECT_EQ(simulation.collect_results()["Sign"].count(),2);
  EXPECT_DOUBLE_EQ(simulation.fraction_completed(),1.);
}

TEST(CtintAtomicOracle, FrequencyTablesIncludeUnmeasuredFrequenciesAndReinitialize) {
  // The independent phase oracle covers table reads beyond the measured prefix,
  // both cached and continuous time, and sequential changes in table dimensions.
  for (bool cached : {false,true,false}) {
    for (unsigned frequencies : {4u,14u,6u}) {
      SCOPED_TRACE(cached);
      SCOPED_TRACE(frequencies);
      auto p=parameters("unused",2.,true);
      p["NMATSUBARA"]=static_cast<int>(frequencies);
      p["NMATSUBARA_MEASUREMENTS"]=2;
      if (cached) p["TAU_DISCRETIZATION_FOR_EXP"]=16;
      c_or_cdagger::initialize_simulation(p);
      for (int sign : {-1,1}) {
        creator op(0,0,.625,frequencies); // Exactly five cells on this tau grid.
        op.compute_exp(frequencies,sign);
        for (unsigned w=0; w<frequencies; ++w) {
          const double phase=(2.*w+1.)*std::acos(-1.)*.625/2.;
          const std::complex<double> expected(std::cos(phase),sign*std::sin(phase));
          EXPECT_NEAR(std::abs(op.exp_iomegat()[w]-expected),0.,1.e-12);
        }
      }
    }
  }
}

TEST(CtintAtomicOracle, MeasuredPrefixWorksWithProductionVertexAndChangedExtents) {
  for (int frequencies : {4,14,6}) {
    SCOPED_TRACE(frequencies);
    auto p=parameters("unused",2.,true);
    p["U"]=4.; p["NMATSUBARA"]=frequencies;
    p["NMATSUBARA_MEASUREMENTS"]=2;
    FrequencyProbe simulation(p,0);
    simulation.measure_one_vertex();
    const auto results=simulation.collect_results();
    EXPECT_EQ(results["Sign"].count(),1);
    for (int flavor=0; flavor<2; ++flavor)
      for (const auto* prefix : {"Wk_real_","Wk_imag_"}) {
        const auto& result=results[std::string(prefix)+std::to_string(flavor)+"_0_0"];
        EXPECT_EQ(result.count(),1);
        const auto values=result.mean<std::vector<double>>();
        ASSERT_EQ(values.size(),2);
        for (double value : values) EXPECT_TRUE(std::isfinite(value));
      }
  }
}

TEST(CtintAtomicOracle, MatrixRebuildPreservesPositiveAndNegativeConfigurationSigns) {
  // Adapted invariant fixture from pinned reconciliation tests/run_config.cpp.
  // This tests a fixed accepted configuration, not a sampled physical average.
  for (double interaction : {-4.,4.}) {
    SCOPED_TRACE(interaction);
    auto p=parameters("unused",2.,true);
    p["U"]=interaction;
    RebuildProbe simulation(p,0);
    simulation.check_rebuild(interaction);
  }
}
