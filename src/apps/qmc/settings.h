/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2001-2006 by Fabien Alet <alet@comp-phys.org>,
*                            Matthias Troyer <troyer@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#ifndef ALPS_APPLICATIONS_QMC_SETTINGS_H
#define ALPS_APPLICATIONS_QMC_SETTINGS_H

#include <alps/expression.h>
#include <alps/parameter.h>

namespace qmc {

// Moved from QMCRun's constructor in src/apps/qmc/qmc.h. Keep the existing alias
// precedence and expression evaluation; constructor validation remains there
// so errors still occur after its lattice/model and measurement initialization.
inline double resolve_inverse_temperature(const alps::Parameters& p)
{
  return p.defined("Beta") ? alps::evaluate<double>(p["Beta"],p)
      : (p.defined("beta") ? alps::evaluate<double>(p["beta"],p)
        : (p.defined("BETA") ? alps::evaluate<double>(p["BETA"],p)
           : (p.defined("T") ? 1./alps::evaluate<double>(p["T"],p)
              : (p.defined("TEMPERATURE") ? 1./alps::evaluate<double>(p["TEMPERATURE"],p)
                 : 1./alps::evaluate<double>(p["temperature"],p)))));
}

} // namespace qmc

#endif // ALPS_APPLICATIONS_QMC_SETTINGS_H
