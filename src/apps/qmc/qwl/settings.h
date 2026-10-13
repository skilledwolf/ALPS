/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2004-2005 by Stefan Wessel <wessel@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#ifndef ALPS_QWL_SETTINGS_H
#define ALPS_QWL_SETTINGS_H

#include <alps/parameter.h>
#include <boost/cstdint.hpp>
#include <cmath>
#include <cstddef>

namespace qwl {

struct wang_landau_schedule {
  boost::uint32_t use_zhou_bhatt;
  unsigned int block_sweeps_total;
  double logf;
  double minimum_histogram;
  double flatness_treshold;
};

// Extracted from QWL_SSE_Simulation's constructor in this repository's
// qwl_sse.h. Preserve its numeric conversions, eager defaults and branch order.
// The caller supplies sizes after allocating the order window and operator
// string, retaining their existing validation and construction order.
inline wang_landau_schedule resolve_wang_landau_schedule(
    const alps::Parameters& parms, unsigned int cutoff,
    unsigned int histogram_size, std::size_t number_sites) {
  wang_landau_schedule schedule;
  schedule.use_zhou_bhatt=parms.value_or_default("USE_ZHOU_BHATT_METHOD",1.);
  if (schedule.use_zhou_bhatt) {
    schedule.block_sweeps_total=histogram_size;
    schedule.logf=std::log(double(parms.value_or_default("INITIAL_MODIFICATION_FACTOR",std::exp(1.))));
    schedule.minimum_histogram=1./schedule.logf;
    schedule.flatness_treshold=parms.value_or_default("FLATNESS_TRESHOLD",1E10);
  }
  else {
    schedule.block_sweeps_total=parms.value_or_default("BLOCK_SWEEPS",10000);
    schedule.logf=std::log(double(parms.value_or_default("INITIAL_INCREASE_FACTOR",std::exp(cutoff*std::log((double)number_sites)/schedule.block_sweeps_total))));
    schedule.minimum_histogram=0.;
    schedule.flatness_treshold=parms.value_or_default("FLATNESS_TRESHOLD",0.2);
  }
  return schedule;
}

} // namespace qwl

#endif
