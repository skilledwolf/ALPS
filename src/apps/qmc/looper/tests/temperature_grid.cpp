// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#include <alps/mc/temperature_grid.hpp>
#include <iostream>
#include <limits>

void require(bool condition,char const* message) {
  if (!condition) throw std::runtime_error(message);
}
template<class F> void rejects(F f) {
  try { f(); } catch (std::invalid_argument const&) { return; }
  throw std::runtime_error("Invalid grid accepted");
}
int main() {
  try {
    alps::params p;
    p["NUM_REPLICAS"]=3; p["BETA_MIN"]=1.; p["BETA_MAX"]=9.; p["TEMPERATURE_DISTRIBUTION_TYPE"]=3;
    alps::mc::temperature_grid root(p);
    require(root.values()==std::vector<double>({1,4,9}),"Square-root beta distribution");
    p.erase("BETA_MIN"); p.erase("BETA_MAX"); p.erase("TEMPERATURE_DISTRIBUTION_TYPE");
    p["T_MIN"]=1.; p["T_MAX"]=3.;
    alps::mc::temperature_grid thermal(p);
    require(thermal.values()==std::vector<double>({1./3,.5,1.}),"Uniform temperature distribution");
    p["TEMPERATURE_SET"]=std::vector<double>{3,1,2};
    require(alps::mc::temperature_grid(p).values()==thermal.values(),"Explicit temperature list");
    rejects([] {alps::mc::temperature_grid grid(std::vector<double>{1,1});});
    rejects([] {alps::mc::temperature_grid grid(std::vector<double>{0,1});});
    rejects([] {alps::mc::temperature_grid grid(std::vector<double>{1,std::numeric_limits<double>::infinity()});});
    rejects([&] {thermal.restore({.25,.5,1});});
    rejects([&] {thermal.restore({1./3,1,.5});});
    auto initial=thermal.values();
    require(!thermal.optimize_population({0,.7,.6}) && thermal.values()==initial,"Invalid feedback mutated temperature grid");
    require(thermal.optimize_population({0,.5,1}),"Uniform diffusive population rejected");
    for (size_t i=0;i<3;++i) require(std::abs(initial[i]-thermal[i])<1e-14,"Uniform diffusion should retain uniform temperatures");
    alps::mc::temperature_grid rate(std::vector<double>{1,2,5});
    rate.optimize_rate(std::vector<double>{1,2,5},[](double energy,double beta){return energy*beta;});
    require(std::abs(rate[1]-3)<.04 && rate[0]==1 && rate[2]==5,"Linear energy law should balance at the midpoint");
    auto optimized=rate.values();
    rate.restore(optimized);
    require(rate.values()==optimized,"Temperature checkpoint altered feedback grid");
    return 0;
  } catch (std::exception const& error) { std::cerr<<error.what()<<'\n'; return 1; }
}
