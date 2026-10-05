// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#include "../qmc.h"
#include "../analysis.hpp"
#include <iostream>

void require(bool value,char const* why) { if (!value) throw std::runtime_error(why); }
class measurements : public QMCRun<> {
public:
  measurements(alps::params const& p,bool signed_values):QMCRun<>(p,8,0) {
    is_signed_=signed_values; initialize_site_states(); create_common_observables();
  }
  void update() override {}
  void measure() override {}
  double fraction_completed() const override { return 0.; }
  bool sample(std::vector<uint8_t> const& states,double sign) { return do_common_measurements(sign,states); }
};
int main() {
  try {
    alps::params p;
    p["LATTICE"]="chain lattice"; p["L"]=4; p["MODEL"]="boson Hubbard";
    p["Nmax"]=2; p["t"]=.3; p["U"]=1.; p["mu"]=.5; p["T"]=1.;
    p["MEASURE[Correlations]"]=true;
    p["MEASURE_CORRELATIONS[nn]"]="n:n";
    p["MEASURE_STRUCTURE_FACTOR[nq]"]="n:n";
    p["INITIAL_SITE"]=1;
    measurements sim(p,true);
    require(sim.sample({0,1,2,1},-1.),"Valid sample rejected");
    auto corr=sim.measurement("nn")->result().mean();
    auto builtin=sim.measurement("Density Correlations")->result().mean();
    for (size_t i=0;i<4;++i) {
      require(corr[i]==-double(std::vector<int>{0,1,2,1}[i]),"Custom origin/sign correlation");
      require(builtin[i]==corr[i],"Built-in and custom correlation disagree");
    }
    require(corr[4]==-1.,"Signed numerator lost its paired denominator");
    auto structure=sim.measurement("nq")->result().mean();
    require(std::abs(structure[0]+4.)<1e-12,"Custom q=0 structure factor normalization");
    p["RESTRICT_MEASUREMENTS[N]"]=3;
    measurements restricted(p,true);
    require(!restricted.sample({0,1,2,1},-1.),"Particle restriction ignored");
    for (auto const& name:restricted.result_names())
      require(std::visit([](auto const& value){return value.count();},restricted.collect_results().at(name))==0,
              "Rejected sample leaked into measurements");
    p.erase("RESTRICT_MEASUREMENTS[N]");
    measurements constant(p,false);
    for (int i=0;i<1000;++i) constant.sample({0,1,2,1},1.);
    auto results=constant.collect_results_as<alps::alea::batch_result<double>>(native_mc::batch_names(constant.get_measurements()));
    native_mc::unavailable_results unavailable;
    native_qmc::derive(results,p,4,unavailable);
    require(results.at("Compressibility").mean()[0]==0. && results.at("Compressibility").stderror()[0]==0.,
            "Constant particle number must have zero compressibility and error");
    return 0;
  } catch (std::exception const& error) { std::cerr<<error.what()<<'\n'; return 1; }
}
