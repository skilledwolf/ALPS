// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
// Shape/finite checks are adapted from read_flavor_vectors in ALPS
// reconciliation 90a36c0d4c64ee7308dcb0651264197ffdf36449,
// src/apps/dmft/qmc/green_function.h (Copyright (C) 2005-2009 Emanuel Gull,
// Philipp Werner, Matthias Troyer and Sebastian Fuchs). This private adapter
// retains both historical CT-INT HDF5 layouts and the NGS parameter interface.
#pragma once

#include "green_function.h"
#include <alps/ngs/params.hpp>
#include <alps/hdf5/vector.hpp>
#include <cmath>
#include <complex>
#include <stdexcept>
#include <string>
#include <vector>

namespace ctint_input {
inline void read_bare_green(const alps::params& parameters,
                            matsubara_green_function_t& green) {
  if (green.nsite()!=1)
    throw std::invalid_argument("CT-INT bare Green function requires one site");
  std::vector<std::vector<std::complex<double>>> values(green.nflavor());
  if (parameters["ATOMIC"] | false) {
    const double beta=parameters["BETA"];
    for (auto& flavor : values) {
      flavor.resize(green.nfreq());
      for (std::size_t w=0; w<flavor.size(); ++w)
        flavor[w]={0.,-beta/((2.*w+1.)*std::acos(-1.))};
    }
  } else {
    alps::hdf5::archive archive(parameters["INFILE"].cast<std::string>(),"r");
    const bool framework=parameters["DMFT_FRAMEWORK"] | false;
    if (framework) {
      unsigned int frequencies, sites, flavors;
      archive["/G0/nt"] >> frequencies;
      archive["/G0/ns"] >> sites;
      archive["/G0/nf"] >> flavors;
      if (frequencies!=green.nfreq() || sites!=1 || flavors!=green.nflavor())
        throw std::invalid_argument("CT-INT /G0 metadata does not match the requested dimensions");
    }
    for (unsigned int f=0; f<green.nflavor(); ++f) {
      const auto path=framework ? "/G0/"+std::to_string(f)+"/mean/value"
                                : "/G0_"+std::to_string(f);
      // This archive API exposes the stored physical dimensions: complex
      // vectors have a trailing real/imaginary axis of length two.
      if (!archive.is_data(path) || !archive.is_complex(path)
          || archive.extent(path)!=std::vector<std::size_t>{green.nfreq(),2})
        throw std::invalid_argument("CT-INT "+path+": expected a complex vector of length "
                                    +std::to_string(green.nfreq()));
      archive[path] >> values[f];
    }
  }
  // Validate every flavor before replacing any caller-owned numerical data.
  for (const auto& flavor : values)
    for (const auto& value : flavor)
      if (!std::isfinite(value.real()) || !std::isfinite(value.imag()))
        throw std::invalid_argument("CT-INT bare Green function must contain only finite values");
  for (unsigned int f=0; f<green.nflavor(); ++f)
    for (unsigned int w=0; w<green.nfreq(); ++w)
      green(w,0,0,f)=values[f][w];
}
}
