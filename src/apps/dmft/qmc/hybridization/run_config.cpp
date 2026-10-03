// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#include <alps/cthyb.hpp>
#include "cthyb_schema.hpp"
#include "input.hpp"
#include <climits>
#include <filesystem>
#include <limits>

namespace alps::cthyb {
std::string_view schema() { return cthyb_schema; }
params prepare_parameters(const params& supplied) {
  auto p=resolve_parameters(supplied,schema());
  if(p["BETA"].as<double>()<=0.) throw std::invalid_argument("BETA must be positive");
  const auto positive=[&](bool enabled,const char* key){
    if(enabled && p[key].as<std::int64_t>()<=0)
      throw std::invalid_argument(std::string(key)+" must be positive for the requested measurement");
  };
  positive(p["MEASURE_freq"].as<bool>() || p["MEASURE_legendre"].as<bool>(),"N_MATSUBARA");
  positive(p["MEASURE_legendre"].as<bool>(),"N_LEGENDRE");
  positive(p["MEASURE_nnt"].as<bool>(),"N_nn");
  positive(p["MEASURE_nnw"].as<bool>(),"N_W");
  const bool two=p["MEASURE_g2w"].as<bool>() || p["MEASURE_h2w"].as<bool>();
  positive(two,"N_w2"); positive(two,"N_W");
  if(two && p["N_w2"].as<std::int64_t>()%2)
    throw std::invalid_argument("N_w2 must be even");
  const auto nw2=p["N_w2"].as<std::int64_t>(), nw=p["N_W"].as<std::int64_t>();
  if(nw2+nw>INT_MAX) throw std::invalid_argument("N_w2 + N_W exceeds the solver index range");
  if(two && (nw2>INT_MAX/nw2 || nw>INT_MAX/(nw2*nw2)))
    throw std::invalid_argument("Two-particle measurement dimensions exceed the solver index range");
  if(p["COMPUTE_VERTEX"].as<bool>()) {
    if(!p["MEASURE_freq"].as<bool>() || !two)
      throw std::invalid_argument("COMPUTE_VERTEX requires MEASURE_freq and MEASURE_g2w or MEASURE_h2w");
    if(p["N_MATSUBARA"].as<std::int64_t>()<nw2/2+nw-1)
      throw std::invalid_argument("COMPUTE_VERTEX requires N_MATSUBARA >= N_w2/2 + N_W - 1");
  }
  const auto orbitals=p["N_ORBITALS"].as<std::int64_t>();
  if(orbitals>INT_MAX/orbitals || p["N_TAU"].as<std::int64_t>()+1>INT_MAX/orbitals)
    throw std::invalid_argument("Orbital/time dimensions exceed the solver index range");
  if(p["MEASURE_sector_statistics"].as<bool>() && orbitals>=31)
    throw std::invalid_argument("Sector statistics require N_ORBITALS < 31");
  if(p["THERMALIZATION"].as<std::uint64_t>()>std::numeric_limits<std::uint64_t>::max()-p["SWEEPS"].as<std::uint64_t>())
    throw std::invalid_argument("Total sweep count exceeds the solver counter range");
  return p;
}
void prepare_run(run_configuration& run) {
  run=resolve_run_configuration(run,schema());
  run.parameters=prepare_parameters(run.parameters);
  const auto orbitals=run.parameters["N_ORBITALS"].as<std::size_t>();
  if(!run.input.exists("interaction_matrix")) {
    if(!run.parameters.exists("U")) throw std::invalid_argument("Specify parameters.U or input.interaction_matrix");
    const double u=run.parameters["U"].as<double>(), j=run.parameters["J"].as<double>();
    const double up=run.parameters.value_or("U'",u-2*j);
    if(!std::isfinite(up)) throw std::invalid_argument("Derived U' must be finite");
    if(orbitals%2 && (up!=u || j!=0.))
      throw std::invalid_argument("Hund interactions require an even N_ORBITALS or an explicit interaction matrix");
  }
  const auto& out=run.output;
  auto base=out["base_path"].as<std::string>();
  if(!base.empty() && base.front()!='/') throw std::invalid_argument("output.base_path must be an absolute HDF5 group path");
  if(base=="/run_config" || base.rfind("/run_config/",0)==0 || base=="/parameters" || base.rfind("/parameters/",0)==0)
    throw std::invalid_argument("output.base_path must not overlap configuration metadata");
  if(base=="/") base.clear();
  if(!base.empty() && base.back()=='/') base.pop_back();
  run.output["base_path"]=base;
  if(run.input["delta_layout"].as<std::string>()=="dmft" && run.input["delta_format"].as<std::string>()!="hdf5")
    throw std::invalid_argument("input.delta_layout='dmft' requires input.delta_format='hdf5'");
  const auto output=std::filesystem::weakly_canonical(out["results"].as<std::string>());
  if(std::filesystem::is_directory(output) || !std::filesystem::is_directory(output.parent_path()))
    throw std::invalid_argument("output.results must name a file in an existing directory");
  for(const auto* key:{"delta","interaction_matrix","chemical_potential","retarded_interaction"})
    if(run.input.exists(key)) {
      const auto filename=run.input[key].as<std::string>();
      if(!std::filesystem::is_regular_file(filename)) throw std::invalid_argument(std::string("Missing input.")+key+" file: "+filename);
      if(std::filesystem::weakly_canonical(filename)==output) throw std::invalid_argument("Output must not replace an input file");
    }
  const auto delta=cthyb_input::series(run.parameters,run.input);
  if(std::any_of(delta.begin(),delta.end(),[](double x){return x>0.;}))
    throw std::invalid_argument("Delta(tau) must be nonpositive");
  if(run.input.exists("interaction_matrix"))
    cthyb_input::static_values(run.input,"interaction_matrix","interaction_format","/Umatrix",orbitals*orbitals);
  if(run.input.exists("chemical_potential"))
    cthyb_input::static_values(run.input,"chemical_potential","chemical_potential_format","/MUvector",orbitals);
  if(run.input.exists("retarded_interaction")) {
    const auto k=cthyb_input::series(run.parameters,run.input,true);
    if(k[0]!=0.) throw std::invalid_argument("Retarded interaction K(tau=0) must be zero");
    for(std::size_t i=0;i<k.size();i+=2)
      if(k[i]<0.) throw std::invalid_argument("Retarded interaction K(tau) must be nonnegative");
  }
  if(out["text"].as<bool>()) {
    const auto directory=std::filesystem::path(out["text_directory"].as<std::string>());
    if(!std::filesystem::is_directory(directory))
      throw std::invalid_argument("output.text_directory must be an existing directory");
    for(const auto* filename:{"simulation.dat","observables.dat","orders.dat","Gt.dat","Ft.dat",
        "Gw.dat","Fw.dat","Sw.dat","Gl_conv.dat","Fl_conv.dat","Gtl.dat","Ftl.dat",
        "Gwl.dat","Fwl.dat","Swl.dat","nnt.dat","nnw.dat","sector_statistics.dat",
        "g2w.dat","h2w.dat","gammaw.dat"}) {
      const auto destination=std::filesystem::weakly_canonical(directory/filename);
      if(destination==output) throw std::invalid_argument("Text output must not replace HDF5 results");
      for(const auto& [key,value]:run.input)
        if((key=="delta" || key=="interaction_matrix" || key=="chemical_potential" || key=="retarded_interaction") &&
           destination==std::filesystem::weakly_canonical(value.as<std::string>()))
          throw std::invalid_argument("Text output must not replace an input file");
    }
  }
}
} // namespace alps::cthyb
