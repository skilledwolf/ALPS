// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#include "../qwl_sse.h"
#include "../evaluation.hpp"
#include <boost/filesystem/operations.hpp>
#include <iostream>

void require(bool value,char const* message) { if (!value) throw std::runtime_error(message); }
void finish(QWL_SSE_Simulation& sim) {
  for (size_t i=0;sim.fraction_completed()<1.;++i) {
    require(i<1000000,"QWL did not converge"); sim.update();
  }
}
void same(QWL_SSE_Simulation const& a,QWL_SSE_Simulation const& b) {
  require(a.completed_sweeps()==b.completed_sweeps(),"QWL restart sweep mismatch");
  auto other=b.collect_results();
  for (auto const& entry:a.collect_results()) std::visit([&](auto const& value) {
    auto const& expected=std::get<std::decay_t<decltype(value)>>(other.at(entry.first));
    require(value.count()==expected.count(),"QWL restart result count mismatch");
    if (value.count()) require(value.mean()==expected.mean(),"QWL restart result mean mismatch");
  },entry.second);
}
int main() {
  auto file=boost::filesystem::temp_directory_path()/boost::filesystem::unique_path("qwl-%%%%-%%%%.h5");
  try {
    alps::params p;
    p["LATTICE"]="chain lattice"; p["MODEL"]="spin"; p["L"]=4; p["J"]=1.;
    p["local_S"]=.5; p["CUTOFF"]=12; p["NUMBER_OF_WANG_LANDAU_STEPS"]=3;
    p["SWEEPS"]=4000; p["SEED"]=137;
    for (auto const* rng:{"mt19937","lagged_fibonacci607"}) {
      p["RNG"]=rng;
      QWL_SSE_Simulation full(p,8),stopped(p,8),resumed(p,8);
      finish(full);
      for (int i=0;i<37;++i) stopped.update();
      { alps::hdf5::archive ar(file.string(),"w"); stopped.save(ar); }
      { alps::hdf5::archive ar(file.string()); resumed.load(ar); }
      same(stopped,resumed); finish(resumed); same(full,resumed);
      require(full.get_random()()==resumed.get_random()(),"QWL restart RNG mismatch");
      // Invalid state must fail before changing either RNG or scientific state.
      { alps::hdf5::archive ar(file.string(),"a"); ar["checkpoint/operators"]<<std::vector<std::array<uint32_t,2>>{{9,0}}; }
      bool rejected=false;
      try { alps::hdf5::archive ar(file.string()); full.load(ar); } catch(std::exception const&) { rejected=true; }
      require(rejected,"Malformed QWL operator checkpoint accepted");
      same(full,resumed);
      require(full.get_random()()==resumed.get_random()(),"Rejected QWL load changed RNG");
    }
    // Exact spin-1/2 AF dimer: H has one -3J/4 and three J/4 eigenvalues.
    // With offset J/4, only the singlet contributes to orders n > 0.
    alps::alea::column<double> g(100); g[0]=std::log(4.);
    for (int n=1;n<100;++n) g[n]=-std::lgamma(n+1.);
    for (double t:{.2,1.,4.}) {
      auto r=qwl::evaluate(g,.25,2,0,t);
      double singlet=std::exp(1/t), z=singlet+3;
      double energy=(.25-singlet/z)/2;
      require(std::abs(r.at("Energy Density")-energy)<1e-12,"QWL dimer energy normalization");
      require(std::abs(r.at("Specific Heat per Site")-3*singlet/(2*t*t*z*z))<1e-12,"QWL dimer heat capacity");
      require(std::abs(r.at("Free Energy Density")-(.25-t*std::log(z))/2)<1e-12,"QWL absolute partition normalization");
    }
    alps::alea::column<double> window(1); window[0]=0;
    auto r=qwl::evaluate(window,1.,4,3,2.);
    require(r.at("Energy Density")==-1.25 && !r.count("Free Energy Density"),"QWL window order/normalization");
    boost::filesystem::remove(file);
    return 0;
  } catch(std::exception const& error) {
    boost::filesystem::remove(file); std::cerr<<error.what()<<'\n'; return 1;
  }
}
