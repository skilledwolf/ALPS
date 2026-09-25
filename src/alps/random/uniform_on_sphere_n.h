/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2002-2005 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#ifndef ALPS_UNIFORM_ON_SPHERE_N_H
#define ALPS_UNIFORM_ON_SPHERE_N_H

#include <boost/random/uniform_on_sphere.hpp>
#include <boost/random/uniform_real.hpp>
#include <cmath> // for std::sqrt
#include <vector>
#include <type_traits>

namespace alps {
namespace detail {

// The low-dimensional algorithms share storage and stateless streaming. The
// general-dimensional distribution remains Boost's normal-based sampler.
template<int N, class RealType, class Cont>
class small_sphere_distribution {
public:
  using input_type = RealType;
  using result_type = Cont;
  explicit small_sphere_distribution(int) : real_(N == 1 ? 0 : -1, 1), container_(N) {}
  void reset() {}

  template<class Engine>
  result_type const& operator()(Engine& eng) {
    if constexpr (N == 1) {
      container_[0] = real_(eng) < 0.5 ? RealType(1) : RealType(-1);
    } else {
      RealType v1, v2, s;
      do {
        v1 = real_(eng);
        v2 = real_(eng);
        s = v1 * v1 + v2 * v2;
      } while (s > 1);
      if constexpr (N == 2) {
        const RealType a = 1.0 / std::sqrt(s);
        container_[0] = a * v1;
        container_[1] = a * v2;
      } else {
        const RealType a = 2 * std::sqrt(1 - s);
        container_[0] = a * v1;
        container_[1] = a * v2;
        container_[2] = 2 * s - 1;
      }
    }
    return container_;
  }

  template<class CharT, class Traits>
  friend std::basic_ostream<CharT, Traits>&
  operator<<(std::basic_ostream<CharT, Traits>& os, small_sphere_distribution const&) {
    return os;
  }
  template<class CharT, class Traits>
  friend std::basic_istream<CharT, Traits>&
  operator>>(std::basic_istream<CharT, Traits>& is, small_sphere_distribution& sd) {
    sd.container_.resize(N);
    return is;
  }

private:
  boost::uniform_real<RealType> real_;
  result_type container_;
};
} // namespace detail

template<int N, class RealType = double, class Cont = std::vector<RealType>>
class uniform_on_sphere_n {
  using base_type = std::conditional_t<(N >= 1 && N <= 3),
    detail::small_sphere_distribution<N, RealType, Cont>,
    boost::uniform_on_sphere<RealType, Cont>>;
public:
  using input_type = typename base_type::input_type;
  using result_type = typename base_type::result_type;
  BOOST_STATIC_CONSTANT(int, dim = N);

  uniform_on_sphere_n() : base_(dim) {}
  void reset() { base_.reset(); }
  template<class Engine>
  result_type const& operator()(Engine& eng) { return base_(eng); }

  template<class CharT, class Traits>
  friend std::basic_ostream<CharT, Traits>&
  operator<<(std::basic_ostream<CharT, Traits>& os, uniform_on_sphere_n const& sd) {
    return os << sd.base_;
  }
  template<class CharT, class Traits>
  friend std::basic_istream<CharT, Traits>&
  operator>>(std::basic_istream<CharT, Traits>& is, uniform_on_sphere_n& sd) {
    return is >> sd.base_;
  }
private:
  base_type base_;
};

} // namespace alps
#endif // ALPS_UNIFORM_ON_SPHERE_N_H
