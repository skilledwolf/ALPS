/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2002-2015 by Matthias Troyer <troyer@comp-phys.org>,
*                            Simon Trebst <trebst@comp-phys.org>,
*                            Synge Todo <wistaria@comp-phys.org>
*                            Lukas Gamper <gamperl@gmail.com>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include "parameter_conversion.hpp"
#include <alps/hdf5/archive.hpp>

namespace {
struct hdf5_writer {
  static const char* extension() { return ".h5"; }
  static const char* name() { return "HDF5"; }
  static void write(const std::string& filename, const alps::Parameters& params) {
    alps::hdf5::archive ar(boost::filesystem::path(filename), "w");
    ar["/parameters"] << params;
  }
};
}

int main(int argc, char** argv) {
  return parameter_conversion::run<hdf5_writer>(argc, argv);
}
