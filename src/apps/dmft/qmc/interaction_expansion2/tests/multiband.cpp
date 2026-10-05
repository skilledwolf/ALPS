// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#include "interaction_expansion.hpp"
#include <alps/ctint.hpp>
#include <Eigen/LU>
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#ifdef ALPS_HAVE_MPI
#include <boost/mpi/environment.hpp>
#endif

namespace {
void require(bool value, char const* message) {
  if (!value) throw std::runtime_error(message);
}
using coupling = std::array<std::array<double, 4>, 4>;
constexpr coupling ragged{{{{0., .7, .5, .3}}, {{.7, 0., .4, 0.}},
                           {{.5, .4, 0., 0.}}, {{.3, 0., 0., 0.}}}};
alps::run_configuration configuration(std::string const& matrix, int sweeps = 8) {
  alps::run_configuration run;
  run.parameters["BETA"] = 2.;
  run.parameters["U"] = 0.;
  run.parameters["MU"] = 0.;
  run.parameters["ALPHA"] = -.01;
  run.parameters["FLAVORS"] = 4;
  run.parameters["N"] = 32;
  run.parameters["NMATSUBARA"] = 5;
  run.parameters["SWEEPS"] = sweeps;
  run.parameters["THERMALIZATION"] = 0;
  run.parameters["MEASUREMENT_PERIOD"] = 1;
  run.input["atomic"] = true;
  run.input["interaction_matrix"] = matrix;
  run.output["results"] = "ctint-multiband-contract.h5";
  run.execution["seed"] = 19;
  run.execution["bins"] = 128;
  alps::ctint::prepare_run(run);
  return run;
}
void write_matrix(std::string const& path, coupling const& interaction) {
  std::ofstream file(path);
  for (std::size_t f=0; f<4; ++f)
    for (std::size_t g=0; g<4; ++g)
      if (interaction[f][g]) file << f << ' ' << g << ' ' << interaction[f][g] << '\n';
  require(bool(file), "cannot write CT-INT coupling fixture");
}

struct fixture : InteractionExpansion {
  using InteractionExpansion::InteractionExpansion;
  void warmup_boundary() {
    for (unsigned i=0; i<3; ++i) { update(); measure(); }
    require(collect_results().at("Sign").count()==0, "last warmup update became a physical measurement");
    update(); measure();
    require(collect_results().at("Sign").count()==1, "first post-warmup update was not measured");
  }
  void check_inverse() const {
    for (std::size_t flavor=0; flavor<M.size(); ++flavor) {
      auto const& block=M[flavor];
      const auto n=num_rows(block.matrix());
      require(n==block.creators().size() && n==block.annihilators().size() && n==block.alpha().size(),
              "vertex removal lost operator bookkeeping");
      if (!n) continue;
      Eigen::MatrixXd direct(n,n);
      for (std::size_t i=0; i<n; ++i)
        for (std::size_t j=0; j<n; ++j) {
          // Independent zero-energy atom: G0(t>0)=-1/2, G0(t<0)=+1/2.
          auto const delta=block.creators()[i].t()-block.annihilators()[j].t();
          direct(i,j)=(delta<0. ? .5 : -.5)+(i==j ? block.alpha()[i] : 0.);
        }
      const Eigen::MatrixXd inverse=direct.inverse();
      double maximum=0.;
      for (std::size_t i=0; i<n; ++i)
        for (std::size_t j=0; j<n; ++j)
          maximum=std::max(maximum,std::abs(block.matrix()(i,j)-inverse(i,j)));
      if (!(maximum<2.e-9)) {
        std::cerr << "inverse flavor=" << flavor << " order=" << n << " max_error=" << maximum
                  << " rcond=" << direct.fullPivLu().rcond() << '\n';
        for (std::size_t i=0; i<n; ++i)
          std::cerr << "operator " << i << " tau=" << block.creators()[i].t()
                    << " alpha=" << block.alpha()[i] << '\n';
        throw std::runtime_error("fast inverse differs from independently inverted atomic determinant");
      }
    }
    for (auto const& entry : vertices) {
      const auto check=[&](unsigned flavor, unsigned index) {
        require(index<M[flavor].creators().size() && M[flavor].creators()[index].flavor()==flavor,
                "vertex points to the wrong flavor or operator");
      };
      check(entry.flavor1(),entry.c_dagger_1());
      check(entry.flavor2(),entry.c_dagger_2());
      require(M[entry.flavor1()].creators()[entry.c_dagger_1()].t()==
              M[entry.flavor2()].creators()[entry.c_dagger_2()].t(), "vertex flavor times disagree");
    }
  }
  void proposals() {
    std::array<unsigned,4> sampled{};
    constexpr unsigned count=4096;
    for (unsigned i=0; i<count; ++i) {
      const auto forward=try_add();
      auto const& entry=vertices.back();
      const auto first=entry.flavor1(), second=entry.flavor2();
      require(first<second && ragged[first][second]!=0., "proposal included a diagonal or absent edge");
      ++sampled[first==0 ? second-1 : 3];
      const auto expected=2.*ragged[first][second]*4.*.75*.75;
      require(std::abs(forward-expected)<1.e-12, "proposal ratio omits the uniform interacting-pair probability");
      if (i<32) {
        perform_add();
        check_inverse();
        const auto reverse=try_remove(0);
        require(std::abs(forward*reverse-1.)<1.e-12, "first-vertex proposal and reverse are not reciprocal");
        perform_remove(0);
      } else reject_add();
    }
    // Degrees (3,2,2,1) expose the old first-flavor/row-conditioned proposal:
    // it favors edge (0,3), with probability 1/3 instead of 1/4.
    require(std::all_of(sampled.begin(),sampled.end(),[](auto count) { return count>900 && count<1150; }),
            "interacting pairs are not proposed uniformly");
    for (unsigned i=0; i<24; ++i) {
      const auto forward=try_add();
      perform_add();
      require(std::abs(forward*try_remove(vertices.size()-1)-1.)<2.e-9,
              "nonempty-configuration proposal and reverse are not reciprocal");
      check_inverse();
    }
    reset_perturbation_series();
    check_inverse();
    while (!vertices.empty()) {
      const auto position=vertices.size()/2;
      try_remove(position); perform_remove(position); check_inverse();
    }
    require(std::all_of(M.begin(),M.end(),[](auto const& block) { return num_rows(block.matrix())==0; }),
            "removing all vertices left determinant operators");
  }
  void noninteracting() {
    for (unsigned i=0; i<17; ++i) { update(); measure(); }
    require(vertices.empty(), "zero matrix created interaction vertices");
    auto const result=collect_results();
    auto const& pairs=result.at("n_i n_j");
    require(pairs.size()==16 && pairs.count()==17, "multiband density-pair extent or cadence changed");
    for (std::size_t f=0; f<4; ++f)
      for (std::size_t g=0; g<4; ++g)
        require(std::abs(pairs.mean()(4*f+g)-(f==g ? .5 : .25))<1.e-12,
                "atomic equal-time pair estimator violates n_f^2=n_f");
  }
  void isolated_flavor() {
    for (unsigned i=0; i<128; ++i) { update(); measure(); }
    require(num_rows(M[3].matrix())==0, "isolated flavor received an interaction vertex");
    require(collect_results().at("densities").mean()(3)==.5, "isolated atomic flavor is no longer free");
  }
};

struct atomic_oracle {
  std::array<double,16> energy{}, weight{};
  double partition=0.;
  explicit atomic_oracle(coupling const& interaction, std::array<double,4> const& levels = {}) {
    for (unsigned state=0; state<16; ++state) {
      for (unsigned f=0; f<4; ++f)
        energy[state]+=levels[f]*double((state>>f)&1);
      for (unsigned f=0; f<4; ++f)
        for (unsigned g=f+1; g<4; ++g)
          energy[state]+=interaction[f][g]*(double((state>>f)&1)-.5)*(double((state>>g)&1)-.5);
      weight[state]=std::exp(-2.*energy[state]);
      partition+=weight[state];
    }
  }
  double pair(unsigned f,unsigned g) const {
    double result=0.;
    for (unsigned state=0; state<16; ++state)
      if (((state>>f)&1) && ((state>>g)&1)) result+=weight[state];
    return result/partition;
  }
  std::complex<double> green(unsigned flavor,unsigned frequency) const {
    const std::complex<double> iw(0.,(2.*frequency+1.)*std::acos(-1.)/2.);
    std::complex<double> result{};
    for (unsigned state=0; state<16; ++state) if (!(state&(1u<<flavor))) {
      const auto occupied=state|(1u<<flavor);
      result+=(weight[state]+weight[occupied])/(iw-(energy[occupied]-energy[state]));
    }
    return result/partition;
  }
  std::array<double,3> moments(unsigned flavor) const {
    std::array<double,3> result{};
    for (unsigned state=0; state<16; ++state) if (!(state&(1u<<flavor))) {
      const auto occupied=state|(1u<<flavor);
      const auto spectral=(weight[state]+weight[occupied])/partition;
      const auto transition=energy[occupied]-energy[state];
      result[0]+=spectral;
      result[1]+=spectral*transition;
      result[2]+=spectral*transition*transition;
    }
    return result;
  }
};

struct tail_fixture : GFourierTransformer {
  using GFourierTransformer::GFourierTransformer;
  std::array<double,5> moments(unsigned f) const {
    return {c1_[f][0][0],c2_[f][0][0],c3_[f][0][0],Sc0_[f][0][0],Sc1_[f][0][0]};
  }
};
void tail_moments(alps::run_configuration run) {
  constexpr std::array<double,4> epsilon{.13,-.21,.37,.08};
  constexpr double mu=.27,field=.09;
  std::array<double,4> levels{};
  run.parameters["MU"]=mu;
  run.parameters["H"]=field;
  for (unsigned f=0; f<4; ++f) {
    levels[f]=epsilon[f]-mu-(f%2 ? field : -field);
    run.parameters["EPS_"+std::to_string(f)]=epsilon[f];
    run.parameters["EPSSQ_"+std::to_string(f)]=epsilon[f]*epsilon[f];
  }
  const atomic_oracle exact(ragged,levels);
  std::vector<double> densities(4),pairs(16);
  for (unsigned f=0; f<4; ++f) {
    densities[f]=exact.pair(f,f);
    for (unsigned g=0; g<4; ++g) pairs[4*f+g]=exact.pair(f,g);
  }
  const U_matrix interaction(run.parameters,run.input);
  const tail_fixture transform(run.parameters,interaction,densities,pairs);
  for (unsigned f=0; f<4; ++f) {
    const auto spectral=exact.moments(f);
    const auto measured=transform.moments(f);
    for (unsigned moment=0; moment<3; ++moment)
      require(std::abs(spectral[moment]-measured[moment])<2.e-12,
              "multiband Fourier moment differs from exact atomic spectral transitions");
    require(std::abs(measured[3]-(spectral[1]-levels[f]))<2.e-12,
            "Hartree tail omits the centered interaction shift");
    require(std::abs(measured[4]-(spectral[2]-spectral[1]*spectral[1]))<2.e-12,
            "self-energy tail omits interflavor density covariance");
  }
}

void atomic_sampling(alps::run_configuration run, coupling const& interaction = ragged) {
  run.parameters["SWEEPS"]=65536;
  run.parameters["THERMALIZATION"]=1024;
  const atomic_oracle exact(interaction);
  fixture simulation(run,0);
  simulation.run([] { return false; });
  auto const result=simulation.collect_results();
  auto const& density=result.at("densities");
  auto const& correlation=result.at("n_i n_j");
  require(density.count()==65536 && correlation.count()==density.count(), "finite atom lost signed samples or included warmup");
  if (interaction[0][3]<0.)
    require(std::abs(result.at("Sign").mean()(0))<.99,
            "mixed-sign atom did not exercise signed normalization");
  for (unsigned f=0; f<4; ++f) {
    require(std::abs(density.mean()(f)-exact.pair(f,f))<8.*density.stderror()(f)+.004,
            "finite atom density differs from occupation-state enumeration");
    for (unsigned g=0; g<4; ++g) {
      const auto index=4*f+g;
      require(std::abs(correlation.mean()(index)-exact.pair(f,g))<8.*correlation.stderror()(index)+.004,
              "finite atom density correlation differs from occupation-state enumeration");
    }
    const auto suffix=std::to_string(f)+"_0_0";
    auto const& real=result.at("Wk_real_"+suffix);
    auto const& imag=result.at("Wk_imag_"+suffix);
    for (unsigned frequency=0; frequency<5; ++frequency) {
      const auto iw=std::complex<double>(0.,(2.*frequency+1.)*std::acos(-1.)/2.);
      const auto bare=1./iw;
      const auto measured=bare-bare*bare*std::complex<double>(real.mean()(frequency),imag.mean()(frequency))/2.;
      const auto error=std::norm(bare)*std::hypot(real.stderror()(frequency),imag.stderror()(frequency))/2.;
      require(std::abs(measured-exact.green(f,frequency))<8.*error+.004,
              "finite atom Green function differs from independent Lehmann sum");
    }
  }
}
}

int main(int argc,char** argv) {
#ifdef ALPS_HAVE_MPI
  boost::mpi::environment environment(argc,argv);
#endif
  const auto matrix=std::filesystem::absolute("ctint-multiband-matrix.dat").string();
  try {
    write_matrix(matrix,ragged);
    const auto run=configuration(matrix);
    auto warmup=run;
    warmup.parameters["THERMALIZATION"]=3;
    fixture(warmup,0).warmup_boundary();
    auto conditioned=run;
    // Forcing every insertion samples configurations that a Markov chain
    // almost never accepts. A larger auxiliary shift keeps this synthetic
    // inverse oracle well conditioned; physical chains retain ALPHA=-.01.
    conditioned.parameters["ALPHA"]=-.25;
    fixture(conditioned,0).proposals();
    tail_moments(run);
    atomic_sampling(run);
    auto mixed=ragged;
    mixed[0][3]=mixed[3][0]=-.15;
    write_matrix(matrix,mixed);
    atomic_sampling(configuration(matrix),mixed);
    auto isolated=ragged;
    isolated[0][3]=isolated[3][0]=0.;
    write_matrix(matrix,isolated);
    fixture(configuration(matrix),0).isolated_flavor();
    write_matrix(matrix,{});
    for (bool histogram : {false,true}) {
      auto free=configuration(matrix);
      free.parameters["HISTOGRAM_MEASUREMENT"]=histogram;
      fixture(free,0).noninteracting();
    }
    // Reusing the API with a different flavor count must not retain the
    // first simulation's static frequency/time work arrays.
    auto two=configuration(matrix);
    two.parameters["FLAVORS"]=2;
    for (unsigned f=2; f<4; ++f) {
      two.parameters.erase("EPS_"+std::to_string(f));
      two.parameters.erase("EPSSQ_"+std::to_string(f));
    }
    alps::ctint::prepare_run(two);
    InteractionExpansion pair(two,0);
    pair.update(); pair.measure();
    require(pair.collect_results().at("densities").size()==2, "two-flavor registry retained multiband extent");
    std::filesystem::remove(matrix);
  } catch (std::exception const& error) {
    std::filesystem::remove(matrix);
    std::cerr << error.what() << '\n';
    return 1;
  }
}
