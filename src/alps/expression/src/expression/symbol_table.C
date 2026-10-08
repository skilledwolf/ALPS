/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2008 by Matthias Troyer <troyer@comp-phys.org>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <alps/expression/symbol_table.h>
#include <iostream>
#include <sstream>

namespace alps {

SymbolTable::SymbolTable(params const& typed)
{
  for (auto const& [key, value] : typed) {
    std::ostringstream text;
    text << value;
    push_back(key, text.str());
  }
}

void SymbolTable::push_back(const entry_type& p, bool allow_overwrite)
{
  if (p.key().empty())
    boost::throw_exception(std::runtime_error("empty key"));
  if (defined(p.key())) {
    if (allow_overwrite)
      map_.find(p.key())->second->value() = p.value();
    else
      boost::throw_exception(std::runtime_error("duplicated parameter: " + p.key()));
  } else {
    list_.push_back(p);
    map_[p.key()] = --list_.end();
  }
}

SymbolTable& SymbolTable::operator<<(const SymbolTable& other)
{
  for (const_iterator it = other.begin(); it != other.end(); ++it)
    (*this) << *it;
  return *this;
}

void SymbolTable::copy_undefined(const SymbolTable& p)
{
  for (const_iterator it=p.begin();it!=p.end();++it)
    if (!defined(it->key()))
      push_back(*it);
}

std::ostream& operator<<(std::ostream& os, const SymbolTable& p)
{
  for (SymbolTable::const_iterator it = p.begin(); it != p.end(); ++it) {
    if (it->value().valid()) {
      std::string s = it->value().c_str();
      os << it->key() << " = ";
      if (s.find(' ') != std::string::npos)
        os << '"' << s << '"';
      else
        os << s;
      os << ";\n";
    }
  }
  return os;
}

} // namespace alps
