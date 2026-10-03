// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#pragma once
#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/vector.hpp>
#include <alps/params.hpp>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

namespace cthyb_input {
inline void check(const std::vector<double>& values, std::size_t size, const std::string& source) {
  if(values.size()!=size) throw std::invalid_argument(source+": data dimensions do not match the configuration");
  if(!std::all_of(values.begin(),values.end(),[](double x){return std::isfinite(x);}))
    throw std::invalid_argument(source+": data must be finite");
}
inline std::vector<double> hdf5_vector(alps::hdf5::archive& ar, const std::string& dataset, std::size_t size) {
  if(ar.extent(dataset)!=std::vector<std::size_t>{size})
    throw std::invalid_argument(dataset+": expected a vector of length "+std::to_string(size));
  std::vector<double> values;
  ar[dataset] >> values;
  check(values,size,dataset);
  return values;
}
inline std::vector<double> text(const std::string& filename, std::size_t rows, std::size_t columns,
                                bool coordinates=false, double beta=0., bool tau=false) {
  std::ifstream stream(filename);
  if(!stream) throw std::invalid_argument("Cannot open input file: "+filename);
  std::vector<double> values(rows*columns);
  for(std::size_t i=0;i<rows;++i){
    if(coordinates){
      double x;
      const double expected=tau?beta*i/(rows-1):double(i);
      if(!(stream>>x)) throw std::invalid_argument(filename+": missing or malformed coordinate data");
      if(!std::isfinite(x) || std::abs(x-expected)>1e-10*std::max(1.,std::abs(expected)))
        throw std::invalid_argument(filename+": expected a uniform "+(tau?std::string("tau"):std::string("index"))+" coordinate grid");
    }
    for(std::size_t j=0;j<columns;++j)
      if(!(stream>>values[i*columns+j])) throw std::invalid_argument(filename+": missing or malformed data");
  }
  stream>>std::ws;
  if(!stream.eof()) throw std::invalid_argument(filename+": unexpected extra data");
  check(values,rows*columns,filename);
  return values;
}
inline std::vector<double> static_values(const alps::params& input, const std::string& key,
                                        const std::string& format_key, const std::string& dataset,
                                        std::size_t size) {
  const auto filename=input[key].as<std::string>();
  if(input[format_key].as<std::string>()=="text") return text(filename,1,size);
  alps::hdf5::archive ar(filename,"r");
  return hdf5_vector(ar,dataset,size);
}
inline std::vector<double> series(const alps::params& parameters, const alps::params& input, bool retarded=false) {
  const std::string key=retarded?"retarded_interaction":"delta";
  const auto filename=input[key].as<std::string>();
  const std::size_t rows=parameters["N_TAU"].as<std::size_t>()+1;
  const std::size_t columns=retarded?2:parameters["N_ORBITALS"].as<std::size_t>();
  const std::string format_key=retarded?"retarded_interaction_format":"delta_format";
  if(input[format_key].as<std::string>()=="text")
    return text(filename,rows,columns,true,parameters["BETA"].as<double>(),
                input[retarded?"retarded_interaction_coordinate":"delta_coordinate"].as<std::string>()=="tau");
  alps::hdf5::archive ar(filename,"r");
  const bool dmft=!retarded && input["delta_layout"].as<std::string>()=="dmft";
  if(dmft){
    unsigned nt,ns,nf;
    ar["/Delta/nt"]>>nt; ar["/Delta/ns"]>>ns; ar["/Delta/nf"]>>nf;
    if(nt!=rows || ns!=1 || nf!=columns) throw std::invalid_argument(filename+": incompatible DMFT Delta dimensions");
  }
  std::vector<double> values(rows*columns);
  for(std::size_t j=0;j<columns;++j){
    const std::string dataset=retarded?(j==0?"/Ret_int_K":"/Ret_int_Kp"):
      (dmft?"/Delta/"+std::to_string(j)+"/mean/value":"/Delta_"+std::to_string(j));
    const auto column=hdf5_vector(ar,dataset,rows);
    for(std::size_t i=0;i<rows;++i) values[i*columns+j]=column[i];
  }
  return values;
}
} // namespace cthyb_input
