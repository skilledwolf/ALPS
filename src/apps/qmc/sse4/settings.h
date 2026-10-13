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
#include <boost/throw_exception.hpp>
#include <stdexcept>
#include <string>

namespace dirloop_sse {

// A sweep or thermalization count, read with StringValue's unsigned integer
// conversion. That conversion wraps a negative value around to a huge count,
// which would never be reached, so negative values are rejected.
inline boost::uint64_t resolve_count(const alps::StringValue& value, const std::string& name)
{
    const boost::uint64_t count = boost::uint64_t(value);
    if (count != 0 && value[0] == '-')
        boost::throw_exception(std::out_of_range("negative " + name + " is illegal"));
    return count;
}

// Extracted from SSE_run's constructor in this directory's sse.h. Keep the
// alias order and StringValue's unsigned integer conversion, including the
// legacy acceptance of zero, separate from the lattice.
inline boost::uint64_t resolve_sweeps(const alps::Parameters& params)
{
    return params.defined("SWEEPS") ?
        resolve_count(params["SWEEPS"], "SWEEPS") : (params.defined("MCS") ?
            resolve_count(params["MCS"], "MCS") : resolve_count(params["Steps"], "Steps"));
}

inline boost::uint64_t resolve_thermalization(const alps::Parameters& params,
                                             boost::uint64_t sweeps)
{
    return params.defined("THERMALIZATION") ?
        resolve_count(params["THERMALIZATION"], "THERMALIZATION") : (params.defined("thermalization") ?
            resolve_count(params["thermalization"], "thermalization") : sweeps / 10);
}

} // namespace dirloop_sse

#endif
