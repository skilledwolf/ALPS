/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2006 by Matthias Troyer <troyer@comp-phys.org>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <alps/parameter/parameter.h>
#include <iostream>

namespace alps {

std::ostream& operator<<(std::ostream& os, const alps::Parameter& p) {
  if (p.value().valid()) {
    std::string s = p.value().c_str();
    os << p.key() << " = ";
    if (s.find(' ') != std::string::npos)
      os << '"' << s << '"';
    else
      os << s;
    os << ";";
  }
  return os;
}

} // namespace alps
