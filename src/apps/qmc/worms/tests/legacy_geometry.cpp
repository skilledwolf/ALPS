// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
// Adapted from reconciliation 98a329a90 geometry tests to the legacy WRun API.
#include "WRun.h"
#include <complex>
#include <boost/numeric/ublas/matrix.hpp>
extern "C" void zheev_(const char*,const char*,const int*,std::complex<double>*,const int*,double*,std::complex<double>*,const int*,double*,int*);
#include <gtest/gtest.h>
#include <complex>
namespace {
alps::Parameters configuration(bool chains=false) {
    alps::Parameters p;
    p["MODEL"]="hardcore boson"; p["LATTICE"]=chains ? "square lattice" : "chain lattice";
    p["L"]=4; if(chains) p["W"]=3;
    p["T"]=.5; p["t"]=1.; p["mu"]=.1; p["V"]=.2;
    p["SEED"]=173; p["DISORDERSEED"]=173; p["RNG"]="mt19937";
    p["THERMALIZATION"]=2000; p["SWEEPS"]=250000;
    p["LATTICE_LIBRARY"]=std::string(ALPS_WORM_LATTICES);
    p["MODEL_LIBRARY"]=std::string(ALPS_WORM_MODELS);
    p["CHAIN_KAPPA"]=chains;
    return p;
}
double ring_free_energy(double flux) {
    // Four hard-core sites, H=-sum (e^(i phi/4)b+_i b_(i+1)+h.c.)
    // +.2 sum n_i n_(i+1)-.1 sum n_i; periodic bonds counted once.
    boost::numeric::ublas::matrix<std::complex<double>,boost::numeric::ublas::column_major> h(16,16,std::complex<double>{});
    const auto phase=std::polar(1.,flux/4);
    for(unsigned state=0;state<16;++state) {
        for(unsigned site=0;site<4;++site) {
            unsigned next=(site+1)%4;
            const bool a=state&(1u<<site), b=state&(1u<<next);
            h(state,state)+=.2*a*b-.1*a;
            if (!a && b) {
                const auto moved=state^(1u<<site)^(1u<<next);
                h(moved,state)-=phase; h(state,moved)-=std::conj(phase);
            }
        }
    }
    const int n=16,lwork=64;int info=0;
    std::vector<double> energies(n),rwork(3*n);
    std::vector<std::complex<double>> work(lwork);
    zheev_("N","U",&n,&h.data()[0],&n,energies.data(),work.data(),&lwork,rwork.data(),&info);
    if(info)throw std::runtime_error("independent ring diagonalization failed");
    const double minimum=energies.front();double partition=0;
    for(double energy:energies)partition+=std::exp(-(energy-minimum)/.5);
    return minimum-.5*std::log(partition);
}
}
TEST(WormGeometry, ClosedWindingMatchesIndependentTwistedRingPartition) {
    auto p=configuration(); p["USE_1D_STIFFNESS"]=true;
    WRun run(alps::ProcessList(1),p,0); run.start_worker();
    for(unsigned i=0;i<252000;++i) run.run();
    const auto& w=run.get_measurements().get<alps::RealVectorObservable>("Winding number^2");
    const auto& histogram=run.get_measurements().get<alps::RealVectorObservable>("Winding number histogram");
    EXPECT_EQ(w.count(),histogram.count());
    const double delta=.002;
    const double exact=(ring_free_energy(delta)-2*ring_free_energy(0)+ring_free_energy(-delta))/(delta*delta);
    EXPECT_GT(exact,.05); EXPECT_NEAR(.5*w.mean()[0],exact,.02);
}
TEST(WormGeometry, ChainSecondMomentsSquareEachMeasuredDensity) {
    auto p=configuration(true); p["THERMALIZATION"]=0;
    WRun run(alps::ProcessList(1),p,0); run.start_worker();
    std::valarray<double> last(0.,3), last_squared(0.,3);
    for(unsigned n=1;n<=200;++n) {
        run.run();
        const auto& density=run.get_measurements().get<alps::RealVectorObservable>("Chain density");
        const auto& second=run.get_measurements().get<alps::RealVectorObservable>("Chain density^2");
        ASSERT_EQ(density.count(),n); ASSERT_EQ(second.count(),n);
        auto sum=density.mean()*double(n); auto squared=second.mean()*double(n);
        ASSERT_EQ(sum.size(),3); ASSERT_EQ(squared.size(),3);
        for(unsigned c=0;c<3;++c) {
            const double value=sum[c]-last[c];
            EXPECT_GE(value,-1.e-12); EXPECT_LE(value,1.+1.e-12);
            EXPECT_NEAR(squared[c]-last_squared[c],value*value,1.e-10);
        }
        last=sum;last_squared=squared;
    }
}
TEST(WormGeometry, OneDimensionalHistogramRejectsHigherDimensions) {
    auto p=configuration(true); p["USE_1D_STIFFNESS"]=true;
    EXPECT_THROW(WRun(alps::ProcessList(1),p,0),std::invalid_argument);
}
