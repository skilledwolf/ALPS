// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
// Extends the HDF5 vector validation cases in ALPS reconciliation
// 90a36c0d4c64ee7308dcb0651264197ffdf36449,
// src/apps/dmft/qmc/interaction_expansion2/tests/run_config.cpp, for the
// retained NGS constructor and both historical input layouts.
#include <gtest/gtest.h>
#include <alps/testing/temporary_directory.hpp>
#include "interaction_expansion2/interaction_expansion.hpp"
#include "interaction_expansion2/bare_green_input.hpp"
#include <limits>

void compute_greens_functions(const alps::results_type<HubbardInteractionExpansion>::type&,
                             const alps::params&, const std::string&);

namespace {
using complex=std::complex<double>;
alps::params parameters(const std::string& file, bool framework) {
  alps::params p;
  p["INFILE"]=file; p["DMFT_FRAMEWORK"]=framework; p["ATOMIC"]=false;
  p["BETA"]=2.; p["U"]=0.; p["MU"]=0.; p["ALPHA"]=-.01;
  p["FLAVORS"]=2; p["SITES"]=1; p["NMATSUBARA"]=4; p["N"]=4;
  p["SWEEPS"]=4; p["THERMALIZATION"]=0; p["MEASUREMENT_PERIOD"]=1;
  p["MAX_TIME"]=60; p["SEED"]=42;
  p["EPS_0"]=0.; p["EPS_1"]=0.; p["EPSSQ_0"]=0.; p["EPSSQ_1"]=0.;
  return p;
}
void write_input(const std::string& file, bool framework,
                 const std::vector<complex>& second) {
  alps::hdf5::archive archive(file,"w");
  if (framework) {
    archive["/G0/nt"] << 4;
    archive["/G0/ns"] << 1;
    archive["/G0/nf"] << 2;
  }
  archive[framework ? "/G0/0/mean/value" : "/G0_0"]
    << std::vector<complex>(4,complex(.125,-.5));
  archive[framework ? "/G0/1/mean/value" : "/G0_1"] << second;
}
class InspectableSimulation : public HubbardInteractionExpansion {
public:
  using HubbardInteractionExpansion::HubbardInteractionExpansion;
  complex frequency(unsigned f, unsigned w) const {
    return bare_green_matsubara(w,0,0,f);
  }
};
}

TEST(CtintBareGreenInput, RejectsWrongLengthBeforeCopyingEitherFlavor) {
  alps::testing::TemporaryDirectory directory;
  const auto file=(directory.path()/"g0.h5").string();
  for (bool framework : {false,true}) {
    SCOPED_TRACE(framework);
    for (std::size_t size : {0,3,5}) {
      SCOPED_TRACE(size);
      write_input(file,framework,std::vector<complex>(size,complex(0.,-.5)));
      const auto p=parameters(file,framework);
      matsubara_green_function_t green(4,1,2);
      for (unsigned f=0; f<2; ++f)
        for (unsigned w=0; w<4; ++w) green(w,0,0,f)={7.,8.};
      EXPECT_THROW(ctint_input::read_bare_green(p,green),std::invalid_argument);
      for (unsigned f=0; f<2; ++f)
        for (unsigned w=0; w<4; ++w) EXPECT_EQ(green(w,0,0,f),complex(7.,8.));
      EXPECT_THROW(HubbardInteractionExpansion(p,0),std::invalid_argument);
    }
  }
}

TEST(CtintBareGreenInput, RejectsNonfiniteAndNoncomplexData) {
  alps::testing::TemporaryDirectory directory;
  const auto file=(directory.path()/"g0.h5").string();
  for (bool framework : {false,true}) {
    SCOPED_TRACE(framework);
    for (complex invalid : {complex(std::numeric_limits<double>::quiet_NaN(),0.),
                            complex(0.,std::numeric_limits<double>::infinity())}) {
      std::vector<complex> values(4,complex(0.,-.5));
      values[3]=invalid;
      write_input(file,framework,values);
      EXPECT_THROW(HubbardInteractionExpansion(parameters(file,framework),0),std::invalid_argument);
    }
    write_input(file,framework,std::vector<complex>(4,complex(0.,-.5)));
    {
      alps::hdf5::archive archive(file,"a");
      archive[framework ? "/G0/1/mean/value" : "/G0_1"] << std::vector<double>(4,-.5);
    }
    EXPECT_THROW(HubbardInteractionExpansion(parameters(file,framework),0),std::invalid_argument);
  }
}

TEST(CtintBareGreenInput, RejectsInconsistentFrameworkMetadata) {
  alps::testing::TemporaryDirectory directory;
  const auto file=(directory.path()/"g0.h5").string();
  for (const auto* field : {"nt","ns","nf"}) {
    SCOPED_TRACE(field);
    write_input(file,true,std::vector<complex>(4,complex(0.,-.5)));
    {
      alps::hdf5::archive archive(file,"a");
      archive[std::string("/G0/")+field] << 3;
    }
    EXPECT_THROW(HubbardInteractionExpansion(parameters(file,true),0),std::invalid_argument);
  }
}

TEST(CtintBareGreenInput, ExplicitAtomicFalseUsesTheFile) {
  alps::testing::TemporaryDirectory directory;
  const auto file=(directory.path()/"g0.h5").string();
  const std::vector<complex> second{{.25,-.5},{.125,-.3},{.0625,-.2},{.03125,-.1}};
  for (bool framework : {false,true}) {
    SCOPED_TRACE(framework);
    write_input(file,framework,second);
    InspectableSimulation simulation(parameters(file,framework),0);
    for (unsigned w=0; w<4; ++w) {
      EXPECT_EQ(simulation.frequency(0,w),complex(.125,-.5));
      EXPECT_EQ(simulation.frequency(1,w),second[w]);
    }
  }
}

TEST(CtintBareGreenInput, AtomicPublicationDoesNotReadAnUnrelatedFile) {
  alps::testing::TemporaryDirectory directory;
  const auto file=(directory.path()/"unrelated.h5").string();
  // Existing input is deliberately malformed; ATOMIC=true must ignore it in
  // both initialization and the production Green-function postprocessor.
  write_input(file,false,std::vector<complex>(1,complex(0.,-.5)));
  auto p=parameters(file,false);
  p["ATOMIC"]=true;
  InspectableSimulation simulation(p,0);
  for (unsigned step=0; step<4; ++step) { simulation.update(); simulation.measure(); }
  const auto previous=std::filesystem::current_path();
  struct RestoreDirectory {
    std::filesystem::path previous;
    ~RestoreDirectory() { std::error_code error; std::filesystem::current_path(previous,error); }
  } restore{previous};
  std::filesystem::current_path(directory.path());
  const auto output=(directory.path()/"result.h5").string();
  compute_greens_functions(simulation.collect_results(),p,output);
  alps::hdf5::archive archive(output,"r");
  for (unsigned f=0; f<2; ++f) {
    std::vector<complex> values;
    archive["/G_omega/"+std::to_string(f)+"/mean/value"] >> values;
    ASSERT_EQ(values.size(),4);
    for (unsigned w=0; w<4; ++w)
      EXPECT_NEAR(std::abs(values[w]-complex(0.,-2./((2.*w+1.)*std::acos(-1.)))),0.,1.e-12);
  }
}
