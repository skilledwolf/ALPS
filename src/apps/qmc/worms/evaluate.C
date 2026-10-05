// Copyright (C) 2002-2003 Matthias Troyer; 2026 ALPS Collaboration.
// SPDX-License-Identifier: MIT
#include "../analysis.hpp"
#include <iostream>

int main(int argc,char** argv) {
  try {
    if (argc<2) throw std::invalid_argument("Usage: worm_evaluate results.h5 [...]");
    for (int i=1;i<argc;++i) {
      std::string filename=argv[i];
      if (filename=="--help") { std::cout<<"worm_evaluate results.h5 [...]\n"; return 0; }
      if (filename.empty() || filename[0]=='-') throw std::invalid_argument("Unknown option: "+filename);
      alps::hdf5::update_archive(filename,[&](auto& ar) {
        alps::params parameters;
        uint64_t sites;
        ar["/parameters"] >> parameters;
        ar["/simulation/number_of_sites"] >> sites;
        std::map<std::string,alps::alea::batch_result<double>> results;
        alps::alea::hdf5_serializer codec(ar,"/simulation/results");
        for (auto const* name:{"Centered Density Moments","Winding number histogram"})
          if (ar.is_group("/simulation/results/"+ar.encode_segment(name)))
            alps::alea::deserialize(codec,ar.encode_segment(name),results[name]);
        native_mc::unavailable_results unavailable;
        native_qmc::derive(results,parameters,sites,unavailable);
        for (auto const* name:{"Compressibility","Superfluid stiffness (1D estimator)"}) {
          if (auto found=results.find(name);found!=results.end()) {
            alps::alea::serialize(codec,ar.encode_segment(name),found->second);
            std::cout<<name<<": "<<found->second.mean()[0]<<" +/- "<<found->second.stderror()[0]<<'\n';
          } else if (unavailable.count(name))
            ar["/simulation/unavailable/"+ar.encode_segment(name)] << unavailable.at(name);
        }
      });
    }
    return 0;
  } catch (std::exception const& error) { std::cerr<<"worm_evaluate: "<<error.what()<<'\n'; return 1; }
}
