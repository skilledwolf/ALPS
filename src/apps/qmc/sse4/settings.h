/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2009-2010 by Sergei Isakov <isakov@itp.phys.ethz.ch>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#ifndef ALPS_APPLICATIONS_DIRLOOP_SSE_SETTINGS_H
#define ALPS_APPLICATIONS_DIRLOOP_SSE_SETTINGS_H

#include <alps/parameter.h>
#include <boost/cstdint.hpp>

namespace dirloop_sse {

// Extracted from SSE_run's constructor in this directory's sse.h. Keep the
// alias order and StringValue's unsigned integer conversion, including the
// legacy acceptance of zero and negative values, separate from the lattice.
inline boost::uint64_t resolve_sweeps(const alps::Parameters& params)
{
    return params.defined("SWEEPS") ?
        boost::uint64_t(params["SWEEPS"]) : (params.defined("MCS") ?
            boost::uint64_t(params["MCS"]) : boost::uint64_t(params["Steps"]));
}

inline boost::uint64_t resolve_thermalization(const alps::Parameters& params,
                                             boost::uint64_t sweeps)
{
    return params.defined("THERMALIZATION") ?
        boost::uint64_t(params["THERMALIZATION"]) : (params.defined("thermalization") ?
            boost::uint64_t(params["thermalization"]) : sweeps / 10);
}

} // namespace dirloop_sse

#endif
