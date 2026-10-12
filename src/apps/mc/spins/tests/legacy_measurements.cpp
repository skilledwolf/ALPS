// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
// Reconciliation 93c26bace/cb58b2745 regressions adapted to original SpinSim.
#include "factory.h"
#include "spinsim.h"
#include <alps/testing/temporary_directory.hpp>
#include <gtest/gtest.h>
#include <fstream>
namespace {
alps::Parameters parameters() {
    alps::Parameters p;
    p["MODEL"]="XY";p["LATTICE"]="dimer";p["UPDATE"]="";
    p["LATTICE_LIBRARY"]=std::string(ALPS_SPINMC_LATTICES);
    p["T"]=1.2;p["J"]=.6;p["S"]=1.;p["h"]=0.;
    p["SWEEPS"]=100000;p["THERMALIZATION"]=0;p["SEED"]=173;p["RNG"]="mt19937";
    return p;
}
template<int N, class Mat=MIdMatrix<double,N>>
struct probe : SpinSim<ONMoment<N>,Mat> {
    explicit probe(const alps::Parameters& p):SpinSim<ONMoment<N>,Mat>(alps::ProcessList(1),p,0) {}
    bool cluster() const {return this->cluster_updates_;}
};
template<int N> void sign_policy() {
    auto p=parameters();p["LATTICE"]="triangle";
    probe<N> ferro(p); EXPECT_TRUE(ferro.cluster());
    p["J"]=-.6;probe<N> antiferro(p);EXPECT_FALSE(antiferro.cluster());
    p["UPDATE"]="cluster";EXPECT_THROW(probe<N>{p},std::runtime_error);
}
template<int N> void staggered_identity() {
    for(double j:{-.6,.6}) {
        auto p=parameters();p["J"]=j;p["UPDATE"]="cluster";
        probe<N> run(p);run.start_worker();run.run();run.run();
        const auto& obs=run.get_measurements();
        const auto moment=alps::RealObsevaluator(obs["|Magnetization|"]);
        const auto staggered=alps::RealObsevaluator(obs["Staggered Magnetization^2"]);
        ASSERT_EQ(moment.count(),1);ASSERT_EQ(staggered.count(),1);
        EXPECT_NEAR(staggered.mean(),1-moment.mean()*moment.mean(),1.e-13);
    }
}
}
TEST(SpinMCLegacyMeasurements, ScalarSignWorksInEvenAndOddDimensions) {
    sign_policy<2>();sign_policy<3>();sign_policy<4>();
}
TEST(SpinMCLegacyMeasurements, FirstStaggeredMomentMatchesExactDimerIdentity) {
    staggered_identity<2>();staggered_identity<3>();staggered_identity<4>();
}
TEST(SpinMCLegacyMeasurements, ZeroEdgesAndSelfBondsDoNotConstrainClusters) {
    alps::testing::TemporaryDirectory dir;const auto file=(dir.path()/"graphs.xml").string();
    {std::ofstream out(file);out<<"<LATTICES><GRAPH name=\"triangle\" vertices=\"3\">"
      "<VERTEX id=\"1\"/><VERTEX id=\"2\"/><VERTEX id=\"3\"/>"
      "<EDGE source=\"1\" target=\"2\" type=\"0\"/>"
      "<EDGE source=\"2\" target=\"3\" type=\"0\"/>"
      "<EDGE source=\"3\" target=\"1\" type=\"1\"/></GRAPH>"
      "<GRAPH name=\"loop\" vertices=\"1\"><VERTEX id=\"1\"/>"
      "<EDGE source=\"1\" target=\"1\"/></GRAPH></LATTICES>";}
    auto p=parameters();p["LATTICE_LIBRARY"]=file;p["LATTICE"]="triangle";
    p["J0"]=-.6;p["J1"]=0.;probe<2> open(p);EXPECT_TRUE(open.cluster());
    p["J1"]=-.6;probe<2> closed(p);EXPECT_FALSE(closed.cluster());
    p=parameters();p["LATTICE_LIBRARY"]=file;p["LATTICE"]="loop";p["J"]=-.6;
    p["UPDATE"]="cluster";probe<2> loop(p);EXPECT_TRUE(loop.cluster());
}
TEST(SpinMCLegacyMeasurements, IsolatedSiteMeasuresBothMatrixRepresentations) {
    alps::testing::TemporaryDirectory dir;const auto file=(dir.path()/"site.xml").string();
    {std::ofstream out(file);out<<"<LATTICES><GRAPH name=\"site\" vertices=\"1\"><VERTEX id=\"1\"/></GRAPH></LATTICES>";}
    auto p=parameters();p["LATTICE_LIBRARY"]=file;p["LATTICE"]="site";
    p["UPDATE"]="local";p["D"]=.7;
    probe<2> scalar(p);probe<2,SquareMatrix<double,2>> square(p);
    scalar.start_worker();square.start_worker();
    for(unsigned n=0;n<10;++n){scalar.run();square.run();}
    EXPECT_NEAR(alps::RealObsevaluator(scalar.get_measurements()["Energy"]).mean(),-.7,1.e-13);
    EXPECT_NEAR(alps::RealObsevaluator(square.get_measurements()["Energy"]).mean(),-.7,1.e-13);
}
