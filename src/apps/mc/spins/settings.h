/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1999-2009 by Matthias Troyer <troyer@comp-phys.org>,
*                            Fabian Stoeckli <fabstoec@student.ethz.ch>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

// Run settings of the classical spin simulation that do not depend on the
// lattice. The parameter interpretation was moved unchanged out of the
// AbstractSpinSim constructor (abstractspinsim.h): same reads, same order,
// same conversions and the same expression evaluation. Couplings, spin
// lengths and onsite terms depend on site and bond types and are still
// resolved by the simulation once its lattice exists.

#ifndef ALPS_APPLICATIONS_MC_SPIN_SETTINGS_H_
#define ALPS_APPLICATIONS_MC_SPIN_SETTINGS_H_

#include <alps/expression.h>
#include <alps/parameter.h>
#include <string>
#include "tinyvec.h"

namespace spinmc {

enum class update_request { automatic, local, cluster, unrecognized };

template <int Dim>
struct settings {
  double beta;
  alps::uint64_t sweeps;            // requested measurement sweeps
  TinyVector<double, Dim> h;
  alps::uint64_t thermalization;    // SWEEPS/10 when THERMALIZATION is absent
  bool quantum_convention;
  double g;
  int print_sweeps;
  bool use_error_limit;             // ERROR_VARIABLE and ERROR_LIMIT both given
  update_request update;
  std::string update_text;          // UPDATE as given, empty when absent
  std::string default_coupling;     // J, used for bond types without J<type>
};

// Reads only through the const interface, so it never inserts parameters.
// Expressions in h and g are evaluated as before, with the default
// evaluator; resolution therefore still draws from the disorder generator
// when they call random functions.
template <int Dim>
settings<Dim> resolve_settings(const alps::Parameters& p)
{
  settings<Dim> s;
  s.beta = p.defined("beta") ? double(p["beta"]) : 1./double(p.required_value("T"));
  s.sweeps = static_cast<alps::uint64_t>(p.required_value("SWEEPS"));
  if (p.defined("h")) {
    // TinyVector's parser takes mutable references but only reads copies.
    std::string text = p["h"];
    alps::Parameters context(p);
    s.h = TinyVector<double, Dim>(text, context);
  }
  else
    s.h = TinyVector<double, Dim>(0.0);
  s.thermalization = alps::uint64_t(p.value_or_default("THERMALIZATION", s.sweeps/10));
  s.quantum_convention = p.defined("CONVENTION") && p["CONVENTION"]=="quantum";
  s.g = p.defined("g") ? alps::evaluate<double>(p["g"], p) : 1.;
  s.print_sweeps = int(p.value_or_default("PRINT_SWEEPS", 0));
  s.use_error_limit = p.defined("ERROR_VARIABLE") && p.defined("ERROR_LIMIT");
  s.update_text = p.defined("UPDATE") ? std::string(p["UPDATE"]) : std::string();
  if (s.update_text == "cluster")
    s.update = update_request::cluster;
  else if (s.update_text == "local")
    s.update = update_request::local;
  else if (s.update_text == "")
    s.update = update_request::automatic;
  else
    s.update = update_request::unrecognized;
  s.default_coupling = p.defined("J") ? std::string(p["J"]) : std::string("1.0");
  return s;
}

} // namespace spinmc

#endif
