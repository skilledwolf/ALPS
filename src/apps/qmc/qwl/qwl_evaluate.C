/*****************************************************************************
* Copyright (C) 2004-2006 Stefan Wessel; 2026 ALPS Collaboration
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*****************************************************************************/
#include "evaluation.hpp"
#include <alps/alea/hdf5.hpp>
#include <alps/params.hpp>
#include <alps/plot.h>
#include <filesystem>
#include <iostream>

void evaluate(std::filesystem::path const& file,std::map<std::string,double> const& options) {
  alps::hdf5::archive ar(file.string());
  alps::params parameters;
  ar["/parameters"] >> parameters;
  uint64_t sites;
  ar["/simulation/number_of_sites"] >> sites;
  auto option=[&](std::string const& name,double fallback) {
    auto found=options.find(name);
    return found==options.end() ? parameters.value_or(name,fallback) : found->second;
  };
  double low=option("T_MIN",.1),high=option("T_MAX",10.),step=option("DELTA_T",.1);
  if (!std::isfinite(low) || !std::isfinite(high) || !std::isfinite(step) || low<=0 || high<low || step<=0 ||
      (high-low)/step>1000000 || (high>low && low+step==low))
    throw std::invalid_argument("Require finite 0 < T_MIN <= T_MAX, DELTA_T > 0 and at most 1000001 temperatures");
  auto first=parameters.value_or<unsigned>("EXPANSION_ORDER_MINIMUM",0);

  struct chain { alps::alea::column<double> coefficients,uniform,staggered; double offset; };
  std::vector<chain> chains;
  auto root="/simulation/realizations/0/clones";
  for (auto const& id:ar.list_children(root)) {
    auto path=std::string(root)+"/"+id;
    bool complete;
    ar[path+"/complete"] >> complete;
    if (!complete) throw std::invalid_argument("QWL evaluation requires completed production on every chain");
    alps::alea::hdf5_serializer codec(ar,path+"/results");
    auto mean=[&](std::string const& name) {
      alps::alea::mean_result<double> value;
      alps::alea::deserialize(codec,ar.encode_segment(name),value);
      if (value.count()!=1) throw std::invalid_argument("Expected one final QWL estimate per chain");
      return value.mean().eval();
    };
    chain value;
    value.coefficients=mean("Coefficients");
    auto offset=mean("Offset");
    if (offset.size()!=1) throw std::invalid_argument("Invalid QWL energy offset");
    value.offset=offset[0];
    if (ar.is_group(path+"/results/Uniform Structure Factor Coefficients")) {
      value.uniform=mean("Uniform Structure Factor Coefficients");
      if (ar.is_group(path+"/results/Staggered Structure Factor Coefficients"))
        value.staggered=mean("Staggered Structure Factor Coefficients");
    }
    chains.push_back(std::move(value));
  }
  if (chains.empty()) throw std::invalid_argument("QWL result has no chains");
  std::map<std::string,alps::plot::Set<double>> curves;
  for (size_t i=0;i<=size_t(std::floor((high-low)/step+.5));++i) {
    double temperature=low+i*step;
    std::map<std::string,alps::alea::mean_acc<double>> estimates;
    for (auto const& chain:chains)
      for (auto const& [name,value]:qwl::evaluate(chain.coefficients,chain.offset,sites,first,temperature,chain.uniform,chain.staggered))
        estimates[name]<<alps::alea::make_adapter(value);
    for (auto const& [name,value]:estimates) curves[name]<<temperature<<value.result().mean()[0];
  }
  const std::map<std::string,std::string> suffixes{{"Energy Density","energy"},{"Free Energy Density","free_energy"},
    {"Entropy Density","entropy"},{"Specific Heat per Site","specific_heat"},
    {"Uniform Structure Factor per Site","uniform_structure_factor"},{"Uniform Susceptibility per Site","uniform_susceptibility"},
    {"Staggered Structure Factor per Site","staggered_structure_factor"}};
  auto prefix=file;
  prefix.replace_extension();
  if (prefix.extension()==".out") prefix.replace_extension();
  for (auto const& [name,curve]:curves) {
    alps::plot::Plot<double> plot(name+" versus Temperature",alps::Parameters(parameters));
    plot.set_labels("Temperature",name); plot<<curve;
    alps::oxstream output(prefix.string()+".plot."+suffixes.at(name)+".xml"); output<<plot;
  }
}

int main(int argc,char** argv) {
  try {
    std::map<std::string,double> options;
    std::vector<std::filesystem::path> files;
    for (int i=1;i<argc;++i) {
      std::string arg=argv[i];
      if (arg=="--help") { std::cout<<"qwl_evaluate [--T_MIN value] [--T_MAX value] [--DELTA_T value] results.h5 [...]\n"; return 0; }
      if (arg=="--T_MIN" || arg=="--T_MAX" || arg=="--DELTA_T") {
        if (++i==argc) throw std::invalid_argument("Missing temperature option value");
        size_t used; std::string value=argv[i]; double number=std::stod(value,&used);
        if (used!=value.size()) throw std::invalid_argument("Invalid temperature option value");
        options[arg.substr(2)]=number;
      } else if (arg.empty() || arg.front()=='-') throw std::invalid_argument("Unknown option: "+arg);
      else files.emplace_back(arg);
    }
    if (files.empty()) throw std::invalid_argument("No QWL result file supplied");
    for (auto const& file:files) evaluate(file,options);
    return 0;
  } catch (std::exception const& error) { std::cerr<<"qwl_evaluate: "<<error.what()<<'\n'; return 1; }
}
