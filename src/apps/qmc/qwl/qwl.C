/*****************************************************************************
* Copyright (C) 2004 Stefan Wessel; 2026 ALPS Collaboration
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*****************************************************************************/
#include "qwl_sse.h"
#include "schema.hpp"
#include "../../mc/native_driver.hpp"

int main(int argc,char** argv) {
  return native_mc::main<QWL_SSE_Simulation>(argc,argv,"qwl",qwl_schema,
    [](std::string const&,toml::node const&) -> char const* { return nullptr; },
    [](alps::params& p,alps::run_configuration const&) {
      if (!p.exists("EXPANSION_ORDER_MAXIMUM")) p["EXPANSION_ORDER_MAXIMUM"]=p["CUTOFF"];
      if (!p.exists("START_STORING")) p["START_STORING"]=p["NUMBER_OF_WANG_LANDAU_STEPS"];
    },
    [](alps::run_configuration const& run,auto const& chains,alps::params const&) {
      std::vector<alps::mcbase::results_type> raw;
      for (auto const& chain:chains) raw.push_back(chain->collect_results());
      alps::mcbase::results_type pooled;
      for (auto const& entry:raw.front()) {
        auto const& name=entry.first;
        std::visit([&](auto const& first) {
        using R=std::decay_t<decltype(first)>;
        std::vector<R> values;
        for (auto const& chain:raw) values.push_back(std::get<R>(chain.at(name)));
        pooled.emplace(name,alps::alea::merge(values));
      },entry.second);
      }
      alps::hdf5::save_checkpoint(run.output["results"].as<std::string>(),[&](alps::hdf5::archive& ar) {
        alps::save_results(pooled,run.parameters,ar,"/simulation/results");
        ar["/run_config"] << run;
        ar["/simulation/number_of_sites"] << uint64_t(chains.front()->number_of_sites());
        ar["/simulation/bipartite"] << chains.front()->bipartite();
        for (size_t id=0;id<chains.size();++id) {
          auto path="/simulation/realizations/0/clones/"+std::to_string(id);
          alps::save_results(raw[id],run.parameters,ar,path+"/results");
          ar[path+"/completed_sweeps"] << chains[id]->completed_sweeps();
          ar[path+"/complete"] << (chains[id]->fraction_completed()==1.);
        }
      });
    });
}
