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
#include <cstdint>

namespace worm {

// Moved from WRun's constructor in Winit.C: required_value is converted directly
// to the original int32_t field type, without expression evaluation. The
// constructor retains its negative-value check after validating eta.
inline std::int32_t resolve_thermalization(const alps::Parameters& p)
{
  return p.required_value("THERMALIZATION");
}

} // namespace worm

#endif // ALPS_APPLICATIONS_WORM_SETTINGS_H
