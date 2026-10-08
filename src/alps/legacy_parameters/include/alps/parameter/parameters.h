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

#ifndef ALPS_PARAMETER_PARAMETERS_H
#define ALPS_PARAMETER_PARAMETERS_H

// for MSVC
#if defined(_MSC_VER)
# pragma warning(disable:4251)
#endif

#include "parameter.h"
#include <alps/params.hpp>
#include <boost/throw_exception.hpp>
#include <iostream>
#include <list>
#include <map>
#include <stdexcept>
#include <string>

/// \file parameters.h
/// \brief the text symbols of the expression, lattice and model libraries

namespace alps {

/// \brief a class storing a set of parameters
///
/// the class acts like an associative array but at the same time remembers the order in which elements were added
class ALPS_DECL Parameters
{
public:
  /// the key (parameter name) is a string
  typedef std::string                     key_type;
  /// the parameter value is a String Value, able to store any type in a text representation
  typedef StringValue                     value_type;
  /// the name-value pair is stored as a Parameter
  typedef Parameter                       parameter_type;
  /// the type of container used internally to store the sequential order
  typedef std::list<parameter_type>       list_type;
  /// an integral type to store the number oif elements
  typedef list_type::size_type            size_type;
  /// \brief the iterator type
  ///
  /// iteration goes in the order of insertion into the class, not alphabetically like in a std::map
  typedef list_type::iterator         iterator;
  /// \brief the const iterator type
  ///
  /// iteration goes in the order of insertion into the class, not alphabetically like in a std::map
  typedef list_type::const_iterator   const_iterator;
  /// the type of container used internally to implment the associative array access
  typedef std::map<key_type, iterator>    map_type;

  /// an empty container of parameters
  Parameters() {}

  /// symbols for the lattice and model libraries, one text value per typed parameter
  explicit Parameters(params const& typed);

  /// copy constructor
  Parameters(Parameters const& params) : list_(params.list_), map_() {
    for (iterator itr = list_.begin(); itr != list_.end(); ++itr) map_[itr->key()] = itr;
  }

  /// assignment operator
  Parameters& operator=(Parameters const& rhs) {
    list_ = rhs.list_;
    map_.clear();
    for (iterator itr = list_.begin(); itr != list_.end(); ++itr) map_[itr->key()] = itr;
    return *this;
  }

  /// erase all parameters
  void clear() { list_.clear(); map_.clear(); }

  /// the number of parameters
  size_type size() const { return list_.size(); }

  /// returns true if size == 0
  bool empty() const { return map_.empty();}

  /// does a parameter with the given name exist?
  bool defined(const key_type& k) const { return (map_.find(k) != map_.end());}

  /// accessing parameters by key (name)
  value_type& operator[](const key_type& k) {
    if (defined(k)) {
      return map_.find(k)->second->value();
    } else {
      push_back(k, value_type());
      return list_.rbegin()->value();
    }
  }

  /// accessing parameters by key (name)
  const value_type& operator[](const key_type& k) const {
    if (!defined(k))
      boost::throw_exception(std::runtime_error("parameter " + k + " not defined"));
    return map_.find(k)->second->value();
  }

  /// \brief erase a parameter with a specific key (this takes O(N) time)
  /// \param k the parameter key (name)
  void erase(key_type const& k) {
    map_type::iterator itr = map_.find(k);
    if (itr != map_.end()) {
      list_.erase(itr->second);
      map_.erase(itr);
    } else {
      std::cerr<<"key not found!"<<std::endl;
    }
  }

  /// \brief returns the value or a default
  /// \param k the key (name) of the parameter
  /// \param v the default value
  /// \return if a parameter with the given name \a k exists, its value is returned, otherwise the default v
  value_type value_or_default(const key_type& k, const value_type& v) const {
    return defined(k) ? (*this)[k] : v;
  }

  /// \brief returns the value of a parameter that must be defined
  /// \param k the key (name) of the parameter
  value_type required_value(const key_type& k) const {
    if (!defined(k))
      boost::throw_exception(std::runtime_error("parameter " + k + " not defined"));
    return map_.find(k)->second->value();
  }

  /// an iterator pointing to the beginning of the parameters
  iterator begin() { return list_.begin(); }
  /// a const iterator pointing to the beginning of the parameters
  const_iterator begin() const { return list_.begin(); }
  /// an iterator pointing past the of the parameters
  iterator end() { return list_.end(); }
  /// a const iterator pointing past the of the parameters
  const_iterator end() const { return list_.end(); }

  /// \brief appends a new parameter to the container
  /// \param p the parameter
  /// \param allow_overwrite indicates whether existing parameters may be overwritten
  /// \throw a std::runtime_error if the parameter key is empty or if it exists already and \a allow_overwrite is false
  void push_back(const parameter_type& p, bool allow_overwrite = false);

  /// \brief appends a new parameter to the container
  /// \param k the parameter key (name)
  /// \param v the parameter value
  /// \param allow_overwrite indicates whether existing parameters may be overwritten
  /// \throw a std::runtime_error if the parameter key \a k is empty or if it exists already and \a allow_overwrite is false
  void push_back(const key_type& k, const value_type& v, bool allow_overwrite = false) {
    push_back(Parameter(k, v), allow_overwrite);
  }

  /// \brief set a parameter value, overwriting any existing value
  Parameters& operator<<(const parameter_type& p) {
    (*this)[p.key()] = p.value();
    return *this;
  }

  /// \brief set parameter values, overwriting any existing value
  Parameters& operator<<(const Parameters& params);

  /// \brief set parameter values, without overwriting existing value
  void copy_undefined(const Parameters& p);

private:
  list_type list_;
  map_type map_;
};

/// write parameters in text-form to a std::ostream
ALPS_DECL std::ostream& operator<<(std::ostream& os, const alps::Parameters& p);

} // namespace alps

#endif // ALPS_PARAMETER_PARAMETERS_H
