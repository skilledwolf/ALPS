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

#ifndef ALPS_EXPRESSION_SYMBOL_TABLE_H
#define ALPS_EXPRESSION_SYMBOL_TABLE_H

// for MSVC
#if defined(_MSC_VER)
# pragma warning(disable:4251)
#endif

#include <alps/config.h>
#include <alps/params.hpp>
#include <alps/stringvalue.h>
#include <boost/throw_exception.hpp>
#include <iostream>
#include <list>
#include <map>
#include <stdexcept>
#include <string>

namespace alps {

/// \brief the named text values of the expression, lattice and model libraries
///
/// A value may be an expression over other symbols, such as J1 = 2*J0, which
/// the libraries evaluate lazily; lattices and models add their own symbols.
/// Applications construct the table from typed alps::params. It acts like an
/// associative array but remembers the order in which symbols were added.
class ALPS_DECL SymbolTable
{
public:
  /// the key (symbol name) is a string
  typedef std::string                     key_type;
  /// the value is a StringValue, able to store any type in a text representation
  typedef StringValue                     value_type;

  /// a named value
  class entry_type
  {
  public:
    entry_type() : key_(), value_() {}
    /// Arbitrary types are converted to text by the StringValue constructor
    template<class U>
    entry_type(const key_type& k, const U& v) : key_(k), value_(v) {}
    key_type& key() { return key_; }
    const key_type& key() const { return key_; }
    value_type& value() { return value_; }
    const value_type& value() const { return value_; }
  private:
    key_type key_;
    value_type value_;
  };

  /// the type of container used internally to store the sequential order
  typedef std::list<entry_type>           list_type;
  /// an integral type to store the number of elements
  typedef list_type::size_type            size_type;
  /// iteration goes in the order of insertion, not alphabetically like in a std::map
  typedef list_type::iterator             iterator;
  /// iteration goes in the order of insertion, not alphabetically like in a std::map
  typedef list_type::const_iterator       const_iterator;
  /// the type of container used internally to implement the associative array access
  typedef std::map<key_type, iterator>    map_type;

  /// an empty table
  SymbolTable() {}

  /// one text value per typed parameter
  explicit SymbolTable(params const& typed);

  SymbolTable(SymbolTable const& other) : list_(other.list_), map_() {
    for (iterator itr = list_.begin(); itr != list_.end(); ++itr) map_[itr->key()] = itr;
  }

  SymbolTable& operator=(SymbolTable const& rhs) {
    list_ = rhs.list_;
    map_.clear();
    for (iterator itr = list_.begin(); itr != list_.end(); ++itr) map_[itr->key()] = itr;
    return *this;
  }

  /// erase all symbols
  void clear() { list_.clear(); map_.clear(); }
  /// the number of symbols
  size_type size() const { return list_.size(); }
  /// returns true if size == 0
  bool empty() const { return map_.empty();}
  /// does a symbol with the given name exist?
  bool defined(const key_type& k) const { return (map_.find(k) != map_.end());}

  /// access by name, adding an empty value if the symbol does not exist
  value_type& operator[](const key_type& k) {
    if (defined(k)) {
      return map_.find(k)->second->value();
    } else {
      push_back(k, value_type());
      return list_.rbegin()->value();
    }
  }

  /// access by name
  /// \throw a std::runtime_error if the symbol does not exist
  const value_type& operator[](const key_type& k) const {
    if (!defined(k))
      boost::throw_exception(std::runtime_error("parameter " + k + " not defined"));
    return map_.find(k)->second->value();
  }

  /// \brief erase a symbol (this takes O(N) time)
  void erase(key_type const& k) {
    map_type::iterator itr = map_.find(k);
    if (itr != map_.end()) {
      list_.erase(itr->second);
      map_.erase(itr);
    } else {
      std::cerr<<"key not found!"<<std::endl;
    }
  }

  /// the value of \a k if it is defined, otherwise \a v
  value_type value_or_default(const key_type& k, const value_type& v) const {
    return defined(k) ? (*this)[k] : v;
  }

  /// the value of a symbol that must be defined
  value_type required_value(const key_type& k) const {
    if (!defined(k))
      boost::throw_exception(std::runtime_error("parameter " + k + " not defined"));
    return map_.find(k)->second->value();
  }

  iterator begin() { return list_.begin(); }
  const_iterator begin() const { return list_.begin(); }
  iterator end() { return list_.end(); }
  const_iterator end() const { return list_.end(); }

  /// \brief appends a new symbol
  /// \param allow_overwrite indicates whether an existing symbol may be overwritten
  /// \throw a std::runtime_error if the key is empty or if it exists already and \a allow_overwrite is false
  void push_back(const entry_type& p, bool allow_overwrite = false);

  /// \brief appends a new symbol
  void push_back(const key_type& k, const value_type& v, bool allow_overwrite = false) {
    push_back(entry_type(k, v), allow_overwrite);
  }

  /// set a value, overwriting any existing value
  SymbolTable& operator<<(const entry_type& p) {
    (*this)[p.key()] = p.value();
    return *this;
  }

  /// set values, overwriting any existing values
  SymbolTable& operator<<(const SymbolTable& other);

  /// set values, without overwriting existing values
  void copy_undefined(const SymbolTable& p);

private:
  list_type list_;
  map_type map_;
};

/// write the symbols as text, one "name = value;" per line
ALPS_DECL std::ostream& operator<<(std::ostream& os, const SymbolTable& p);

} // namespace alps

#endif // ALPS_EXPRESSION_SYMBOL_TABLE_H
