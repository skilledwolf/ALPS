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

// Spin lengths, onsite terms and couplings of the classical spin simulation.
// The parameter interpretation was moved unchanged out of the AbstractSpinSim
// constructor (abstractspinsim.h): same reads, same order, same fallbacks and
// the same expression evaluation. The lattice enters only through the site
// types and bond classes it contains, listed in the order the simulation
// visits them, so the terms resolve without building a lattice graph.
//
// The parameter names recorded for each term are new; they report which
// parameter supplied a value and do not change what is read. Also new: J' and
// J'' are rejected for bond types 1 and 2 without J1 or J2 (see below).

#ifndef ALPS_APPLICATIONS_MC_SPIN_COUPLINGS_H_
#define ALPS_APPLICATIONS_MC_SPIN_COUPLINGS_H_

#include <alps/expression.h>
#include <alps/parameter.h>
#include <boost/lexical_cast.hpp>
#include <boost/throw_exception.hpp>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace spinmc {

// Terms of one site type.
template <class MAT>
struct site_terms {
  int type;
  std::string spin_parameter;          // S<type> or S; empty for the default
  double spin;                         // 1, or 1/2 in the quantum convention
  double spinfactor;                   // S, or sqrt(S(S+1)) in the quantum convention
  std::string onsite_parameter;        // D<type> or D; empty when absent
  MAT selfinteraction;                 // D scaled by S^2
};

// Bonds with the same bond type and the same source and target site types.
struct bond_class {
  int type;
  int source_type;
  int target_type;
};

// Terms of one bond class. The simulation renumbers bond types so that a
// class's position in the list is its new bond type.
template <class MAT>
struct bond_terms {
  bond_class bonds;
  std::string coupling_parameter;      // J<type> or J; empty for the default 1.0
  MAT coupling;                        // -J, or J in the quantum convention
  double det;
};

template <class MAT>
struct couplings {
  std::vector<site_terms<MAT> > sites;
  std::vector<bond_terms<MAT> > bonds;
  bool general_case;                   // some site type has an onsite term
  bool ferromagnetic;                  // no coupling has a positive determinant
  bool antiferromagnetic;              // no coupling has a negative determinant
};

// default_coupling is the text used for a bond type without J<type>: the
// value of J, or "1.0" (settings<Dim>::default_coupling).
template <class MAT>
couplings<MAT> resolve_couplings(const alps::Parameters& parms,
                                 const std::vector<int>& site_types,
                                 const std::vector<bond_class>& bond_classes,
                                 bool quantum_convention,
                                 const std::string& default_coupling)
{
  couplings<MAT> result;
  result.general_case = false;
  result.ferromagnetic = true;
  result.antiferromagnetic = true;

  for (std::size_t i = 0; i < site_types.size(); ++i) {
    site_terms<MAT> site;
    site.type = site_types[i];

    // set the values for the spin factors
    std::string spin = "S"+boost::lexical_cast<std::string, int>(site.type);
    if (!parms.defined(spin))
      spin="S";
    if (!parms.defined(spin))
      site.spin = quantum_convention ? 0.5 : 1.;
    else {
      site.spin = alps::evaluate<double>(parms[spin],parms);
      site.spin_parameter = spin;
    }
    site.spinfactor=std::sqrt(site.spin*(site.spin+(quantum_convention ? 1. : 0.)));

    // set the values for the diagonal matrix terms
    std::string diag="D"+boost::lexical_cast<std::string, int>(site.type);
    if (!parms.defined(diag)) {
      diag="D";
    }
    if (!parms.defined(diag))
      site.selfinteraction.setMatrix(0.0);
    else {
      // parse input string to obtain values for the diagonal term
      std::string diag_string = parms[diag];
      site.selfinteraction.setMatrix(diag_string,parms);
      site.onsite_parameter = diag;
      result.general_case=true;
    }
    site.selfinteraction = site.selfinteraction*site.spin*site.spin;

    if (!site.selfinteraction.is_symmetric())
       boost::throw_exception(std::runtime_error("Invalid matrix for D"));
    result.sites.push_back(site);
  }

  for (std::size_t i = 0; i < bond_classes.size(); ++i) {
    bond_terms<MAT> bond;
    bond.bonds = bond_classes[i];
    const std::string name = "J"+boost::lexical_cast<std::string,int>(bond.bonds.type);
    // The model library names the couplings of bond types 1 and 2 J' and J''.
    // Only J<type> and J are read here, so reject a primed coupling that would
    // otherwise be ignored in favour of J.
    const std::string primed = bond.bonds.type == 1 ? "J'" : bond.bonds.type == 2 ? "J''" : "";
    if (!primed.empty() && parms.defined(primed) && !parms.defined(name))
      boost::throw_exception(std::runtime_error(
        primed + " is not read: bond type " + boost::lexical_cast<std::string,int>(bond.bonds.type)
        + " takes " + name + ", falling back to J. Set " + name + ", for example " + name
        + "=\"" + primed + "\", or remove " + primed));
    std::string J_name = parms.value_or_default(name,default_coupling);
    MAT this_J;
    this_J.setMatrix(J_name,parms);
    if (parms.defined(name))
      bond.coupling_parameter = name;
    else if (parms.defined("J"))
      bond.coupling_parameter = "J";
    this_J = this_J * (quantum_convention ? 1. : -1.);
    bond.coupling=this_J;
    double det = this_J.det();
    bond.det = det;
    if (det>0.)
      result.ferromagnetic=false;
    else if (det<0.)
      result.antiferromagnetic=false;
    result.bonds.push_back(bond);
  }
  return result;
}

} // namespace spinmc

#endif
