/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1997-2010 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include "application.hpp"
#include "analysis.hpp"
#include "../native_driver.hpp"
#include "schema.hpp"

int main(int argc,char** argv) {
  auto schema=native_qmc::schema("loop",qmc_common_schema,qmc_application_schema)+alps::mc::replica_exchange_schema;
  return alps::mc::main<looper::application>(argc,argv,"loop",schema.c_str(),
    [](std::string const&)->char const*{return nullptr;},
    [](alps::params& p,alps::run_configuration const& run) {
      if (!p.exists("THERMALIZATION")) p["THERMALIZATION"]=p["SWEEPS"].as<uint64_t>()/10;
      if ((p.value_or("OPTIMIZE_TEMPERATURE",false) || p.value_or("TEMPERATURE_OPTIMIZATION",false)) && run.execution["chains"].as<size_t>()!=1)
        throw std::invalid_argument("Temperature optimization requires one replica ladder (execution.chains = 1)");
    }, [](auto const& run,auto const& chains,auto const&) {
      alps::hdf5::save_checkpoint(run.output["results"].template as<std::string>(),[&](alps::hdf5::archive& ar) {
        for (size_t i=0;i<chains.front()->replicas();++i) {
          std::vector<looper::application::replica_view> views;
          for (auto const& chain:chains) views.push_back({*chain,i});
          std::vector<looper::application::replica_view const*> pointers;
          for (auto const& view:views) pointers.push_back(&view);
          auto replica_run=run;
          if (chains.front()->ladder()) {
            replica_run.parameters["T"]=1/chains.front()->temperatures()[i];
            replica_run.parameters.erase("BETA");
          }
          auto path=chains.front()->ladder() ? "/simulation/replicas/"+std::to_string(i) : "/simulation";
          native_qmc::publish_into(ar,replica_run,pointers,path,looper::derive);
          ar[path+"/parameters"] << replica_run.parameters;
        }
        ar["/parameters"] << run.parameters;
        ar["/run_config"] << run;
      });
    });
}
