/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2012 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <alps/random/uniform_on_sphere_n.h>
#include <alps/utility/vectorio.hpp>
#include <boost/random.hpp>
#include <iostream>
#include <vector>
#include <deque>
#include <limits>
#include <sstream>
#include <stdexcept>

void require(bool condition) {
  if (!condition) throw std::runtime_error("sphere sampling contract failed");
}

template<int N, class Real, class Container = std::vector<Real>>
void check_distribution() {
  boost::mt19937 engine(729);
  alps::uniform_on_sphere_n<N, Real, Container> distribution;
  for (int i = 0; i < 1000; ++i) {
    auto const& sample = distribution(engine);
    require(sample.size() == N);
    Real norm = 0;
    for (auto x : sample) norm += x*x;
    require(std::isfinite(norm) && std::abs(norm - 1) < 20*std::numeric_limits<Real>::epsilon());
  }
  std::stringstream state;
  state << distribution;
  alps::uniform_on_sphere_n<N, Real, Container> restored;
  state >> restored;
  auto saved_engine = engine;
  require(distribution(engine) == restored(saved_engine));
  require(engine == saved_engine);
  distribution.reset();
}

struct scripted_engine {
  using result_type = double;
  static constexpr double min() { return 0; }
  static constexpr double max() { return 1; }
  double operator()() { return values.at(index++); }
  std::vector<double> values;
  std::size_t index = 0;
};

void check_low_dimensions() {
  scripted_engine signs{{.25, .75}};
  alps::uniform_on_sphere_n<1> one;
  require(one(signs)[0] == 1 && one(signs)[0] == -1 && signs.index == 2);
  scripted_engine circle{{.99, .99, .65, .7}};
  alps::uniform_on_sphere_n<2> two;
  auto xy = two(circle);
  require(circle.index == 4 && std::abs(xy[0] - .6) < 1e-14 && std::abs(xy[1] - .8) < 1e-14);
  scripted_engine sphere{{.99, .99, .65, .7}};
  alps::uniform_on_sphere_n<3> three;
  auto xyz = three(sphere);
  require(sphere.index == 4 && std::abs(xyz[0] - .6*std::sqrt(.75)) < 1e-14);
  require(std::abs(xyz[1] - .8*std::sqrt(.75)) < 1e-14 && std::abs(xyz[2] + .5) < 1e-14);
}

int main() {
  check_low_dimensions();
  check_distribution<1, double>();
  check_distribution<2, double>();
  check_distribution<3, double>();
  check_distribution<4, double>();
  check_distribution<5, double>();
  check_distribution<1, float>();
  check_distribution<2, float>();
  check_distribution<3, float>();
  check_distribution<4, float>();
  check_distribution<3, double, std::deque<double>>();
  static const int dim = 3;
  boost::mt19937 eng;
  alps::uniform_on_sphere_n<dim, double, std::vector<double> > dist;
  for (int i = 0; i < 100; ++i) {
    std::vector<double> r = dist(eng);
    std::cout << '(' << alps::write_vector(r, ",", 3) << ')' << std::endl;
  }
  return 0;
}
