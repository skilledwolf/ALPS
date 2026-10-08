/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2009 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#ifndef ALPS_PARAMETER_PARAMETER_H
#define ALPS_PARAMETER_PARAMETER_H

#include <alps/config.h>
#include <alps/stringvalue.h>
#include <iosfwd>
#include <string>

/// \file parameter.h
/// \brief a named text value of the expression, lattice and model libraries

namespace alps {

/// \brief a class to store a single parameter value
///
/// the parameter name (key) is stored as a std::string
/// the parameter value is stored as a StringValue.
class ALPS_DECL Parameter
{
public:
  /// the parameter name (key) is stored as a std::string
  typedef std::string key_type;
  /// the parameter value is stored as a StringValue.
  typedef StringValue value_type;

  /// deault constructor: no name and no value
  Parameter() : key_(), value_() {}
  /// \brief a parameter with a name and value.
  ///
  /// Arbitrary types can be stored. The StringValue constructor will convert
  /// them to a string using boost::lexical_cast
  template<class U>
  Parameter(const key_type& k, const U& v) : key_(k), value_(v) {}

  /// returns the key (parameter name)
  key_type& key() { return key_; }
  /// returns the key (parameter name)
  const key_type& key() const { return key_; }
  /// returns the value
  value_type& value() { return value_; }
  /// returns the value
  const value_type& value() const { return value_; }

private:
  key_type key_;
  value_type value_;
};

/// write parameter in text-form to a std::ostream
ALPS_DECL std::ostream& operator<<(std::ostream& os, const alps::Parameter& p);

} // namespace alps

#endif // ALPS_PARAMETER_PARAMETER_H
