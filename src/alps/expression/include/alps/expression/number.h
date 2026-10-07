/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2010 by Matthias Troyer <troyer@comp-phys.org>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#ifndef ALPS_EXPRESSION_NUMBER_H
#define ALPS_EXPRESSION_NUMBER_H

#include <alps/expression/expression_fwd.h>
#include <alps/expression/evaluate_helper.h>
#include <alps/type_traits/real_type.hpp>
#include <alps/utility/round_trip.hpp>
#include <boost/call_traits.hpp>

namespace alps {
namespace expression {

template<class T>
class Number : public Evaluatable<T> {
public:
  typedef T value_type;
  typedef typename alps::real_type<T>::type real_type;

  Number(typename boost::call_traits<value_type>::param_type x) : val_(x) {}
  value_type value(const Evaluator<T>& =Evaluator<T>(), bool=false) const;
  bool can_evaluate(const Evaluator<T>& =Evaluator<T>(), bool=false) const { return true; }
  void output(std::ostream&) const;
  Evaluatable<T>* clone() const { return new Number<T>(*this); }
private:
  value_type val_;
};

template<class T>
typename Number<T>::value_type Number<T>::value(const Evaluator<T>&, bool) const
{
  return val_;
}

// Expressions are reparsed from their text, so numbers print exactly.
template<class T>
void Number<T>::output(std::ostream& os) const
{
  if (evaluate_helper<T>::imag(val_) == 0)
    os << round_trip_string(evaluate_helper<T>::real(val_));
  else
    os << '(' << round_trip_string(evaluate_helper<T>::real(val_)) << ','
       << round_trip_string(evaluate_helper<T>::imag(val_)) << ')';
}

} // end namespace expression
} // end namespace alps

#endif // ! ALPS_EXPRESSION_IMPL_H
