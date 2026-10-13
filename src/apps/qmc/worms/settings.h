/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2001-2005 by Matthias Troyer <troyer@comp-phys.org>,
*                            Simon Trebst <trebst@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#ifndef ALPS_APPLICATIONS_WORM_SETTINGS_H
#define ALPS_APPLICATIONS_WORM_SETTINGS_H

#include <alps/parameter.h>
#include <boost/throw_exception.hpp>
#include <cmath>
#include <cstdint>
#include <stdexcept>

namespace worm {

// Moved from WRun's constructor in Winit.C: required_value is converted directly
// to the original int32_t field type, without expression evaluation. The
// constructor retains its negative-value check after validating eta.
inline std::int32_t resolve_thermalization(const alps::Parameters& p)
{
  return p.required_value("THERMALIZATION");
}

// The number of sweeps after thermalization, converted to a double without
// expression evaluation, as WRun::work_done always did. A negative count made
// the progress negative and an infinite one kept it at zero, so the run never
// finished; both are rejected, as is a count that is not a number. Zero is
// accepted, as before.
inline double resolve_sweeps(const alps::Parameters& p)
{
  const double sweeps = p.required_value("SWEEPS");
  if (!std::isfinite(sweeps))
    boost::throw_exception(std::invalid_argument("SWEEPS must be a finite number"));
  if (sweeps < 0)
    boost::throw_exception(std::out_of_range("negative SWEEPS is illegal"));
  return sweeps;
}

} // namespace worm

#endif // ALPS_APPLICATIONS_WORM_SETTINGS_H
