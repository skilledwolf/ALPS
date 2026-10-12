// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
// Adapted from reconciliation d0478b0b5 to the retained SSE_run interface.
#include "sse.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
extern "C" void dsyev_(const char*,const char*,const int*,double*,const int*,double*,double*,const int*,int*);
#include <gtest/gtest.h>
namespace {
alps::Parameters configuration(const std::string& family,unsigned skip,bool triangle,bool origin,bool free) {
    alps::Parameters p;
    p["MODEL"]=family;p["LATTICE"]="dimer";p["L"]=2;
    p["SWEEPS"]=200000*skip;p["THERMALIZATION"]=2000;p["SKIP"]=skip;
    p["SEED"]=173;p["RNG"]="mt19937";p["DISORDERSEED"]=173;p["T"]=1.2;
    p["INITIAL_CUTOFF"]=128;p["NUMBER_OF_WORMS_PER_SWEEP"]=7;p["WHICH_LOOP_TYPE"]="heatbath";
    if(family=="spin"){p["local_S"]=1.;p["J"]=.6;p["h"]=.3;p["D"]=.15;}
    else {p["t"]=.45;p["U"]=.8;p["V"]=.2;p["mu"]=.3;if(family=="boson Hubbard")p["Nmax"]=2;}
    p["LATTICE_LIBRARY"]=std::string(ALPS_DIRLOOP_LATTICES);
    p["MODEL_LIBRARY"]=std::string(ALPS_DIRLOOP_MODELS);
    p["MEASURE[Green Function]"]=true;if(origin)p["INITIAL_SITE"]=0;
    if(triangle){p["LATTICE"]="chain lattice";p["L"]=3;p["local_S"]=.5;p["J"]=1.;p["h"]=0.;p["D"]=0.;}
    if(free){for(const char* key:{"t","U","V","mu"})p[key]=0.;p["EPSILON"]=.1;}
    return p;
}
inline std::vector<double> independent_green(const std::string &family,bool triangle=false,bool free=false) {
    const bool spin=family=="spin";const unsigned sites=triangle?3:2;
    const unsigned states=triangle || family=="hardcore boson"?2:3;
    const double magnitude=triangle?.5:1.;
    unsigned dimension=1;for(unsigned site=0;site<sites;++site)dimension*=states;
    namespace ublas=boost::numeric::ublas;
    using matrix=ublas::matrix<double,ublas::column_major>;
    matrix h(dimension,dimension,0.);
    std::vector<matrix> raise(sites,h),charge(sites,h);
    for(unsigned column=0;column<dimension;++column) {
        unsigned remainder=column,stride=1;
        for(unsigned site=0;site<sites;++site,stride*=states) {
            const unsigned level=remainder%states;remainder/=states;
            const double m=double(level)-magnitude;
            charge[site](column,column)=spin?m:double(level);
            if(spin)h(column,column)+=triangle?0:-.3*m+.15*m*m;
            else if(!free)h(column,column)+=.4*level*(double(level)-1)-.3*level;
            if(level+1<states)raise[site](column+stride,column)=spin
                ? std::sqrt(magnitude*(magnitude+1)-m*(m+1)):std::sqrt(double(level+1));
        }
    }
    const std::vector<std::pair<unsigned,unsigned>> edges=triangle
        ? std::vector<std::pair<unsigned,unsigned>>{{0,1},{1,2},{2,0}}
        : std::vector<std::pair<unsigned,unsigned>>{{0,1}};
    for(const auto &edge:edges) {
        const auto first=edge.first,second=edge.second;
        if(spin)h+=(triangle?1.:.6)*ublas::prod(charge[first],charge[second]);
        else if(!free)h+=.2*ublas::prod(charge[first],charge[second]);
        const double transfer=spin?(triangle?.5:.3):free?0:-.45;
        h+=transfer*(ublas::prod(raise[first],ublas::trans(raise[second]))+ublas::prod(ublas::trans(raise[first]),raise[second]));
    }
    const int n=dimension,lwork=4*n;int info=0;
    std::vector<double> energies(n),work(lwork);
    dsyev_("V","U",&n,&h.data()[0],&n,energies.data(),work.data(),&lwork,&info);
    if(info)throw std::runtime_error("Independent Green Hamiltonian trace failed");
    std::vector<double> green(sites*sites,0.);long double partition=0;
    for(unsigned state=0;state<dimension;++state) {
        const long double weight=std::exp(-(energies[state]-energies[0])/1.2);partition+=weight;
        const ublas::vector<double> vector=ublas::column(h,state);
        for(unsigned first=0;first<sites;++first)for(unsigned second=0;second<sites;++second) {
            // Gij=.5<lower_i raise_j + raise_i lower_j>, including the
            // actual truncated onsite commutator (not <n>+.5 at top cutoff).
            const matrix op=.5*(ublas::prod(ublas::trans(raise[first]),raise[second])+ublas::prod(raise[first],ublas::trans(raise[second])));
            green[first*sites+second]+=double(weight)*ublas::inner_prod(vector,ublas::prod(op,vector));
        }
    }
    for(auto &value:green)value/=double(partition);
    return green;
}
inline std::vector<double> output_trace(const alps::Parameters &parameters,
    const std::vector<double> &all_pairs,bool origin) {
    alps::graph_helper<> graph(parameters);const auto sites=graph.num_sites();
    if(origin)return std::vector<double>(all_pairs.begin(),all_pairs.begin()+sites);
    std::vector<double> result(graph.num_distances(),0.);std::vector<unsigned> count(result.size(),0);
    for(unsigned first=0;first<sites;++first)for(unsigned second=0;second<sites;++second) {
        const auto distance=graph.distance(first,second);result.at(distance)+=all_pairs.at(first*sites+second);++count.at(distance);
    }
    for(unsigned distance=0;distance<result.size();++distance)result[distance]/=count[distance];
    return result;
}
void check_green(const std::string& family,unsigned skip=3,bool triangle=false,bool origin=true,bool free=false) {
    auto p=configuration(family,skip,triangle,origin,free);
    SSE_run run(alps::ProcessList(1),p,0);run.start_worker();
    for(unsigned i=0;i<2000+200000*skip;++i)run.run();
    const auto exact=output_trace(p,independent_green(family,triangle,free),origin);
    const alps::RealVectorObsevaluator actual(run.measurements()["Green's Function"]);
    // Signed ratio evaluators retain only complete reblocked bins; verify
    // submissions on the original observable before evaluating the ratio.
    const auto& raw=run.measurements().get<alps::AbstractSimpleObservable<std::valarray<double>>>("Green's Function");
    ASSERT_EQ(raw.count(),200000);ASSERT_EQ(actual.mean().size(),exact.size());
    for(unsigned i=0;i<exact.size();++i){SCOPED_TRACE(i);EXPECT_NEAR(actual.mean()[i],exact[i],.035);}
}
}
TEST(SseLegacyGreen, SignedOpenWormWithSkippedMeasurements){check_green("spin");}
TEST(SseLegacyGreen, EverySweepNormalization){check_green("spin",1);}
TEST(SseLegacyGreen, DistanceClassNormalization){check_green("spin",3,false,false);}
TEST(SseLegacyGreen, HardcoreDimer){check_green("hardcore boson");}
TEST(SseLegacyGreen, FiniteBosonDimer){check_green("boson Hubbard");}
TEST(SseLegacyGreen, FreeTruncatedBosonLadders){check_green("boson Hubbard",3,false,true,true);}
TEST(SseLegacyGreen, FrustratedTriangle){check_green("spin",3,true);}
TEST(SseLegacyGreen, RejectsIncompleteWormWeights){
    for(const char* key:{"WORM_ABORT","RESTRICT_MEASUREMENTS[N]","RESTRICT_MEASUREMENTS[Sz]"}){
        auto p=configuration("spin",1,false,true,false);p[key]=1;
        EXPECT_THROW(SSE_run(alps::ProcessList(1),p,0),std::invalid_argument);
    }
}
