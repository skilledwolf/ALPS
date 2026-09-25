/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2009 by Emanuel Gull <gull@phys.columbia.edu>
 *                              Philipp Werner <werner@itp.phys.ethz.ch>,
 *                              Sebastian Fuchs <fuchs@theorie.physik.uni-goettingen.de>
 *                              Matthias Troyer <troyer@comp-phys.org>
 *
 *
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include "types.hpp"
#include <ostream>

std::ostream& operator<<(std::ostream& os, vertex_array const& vertices) {
  for (auto const& v : vertices) os << v << std::endl;
  return os;
}

std::ostream& operator<<(std::ostream& os, vertex const& v) {
  return os << "z1: " << v.flavor1() << " c1^dag: " << v.c_dagger_1()
            << " c1: " << v.c_1() << " z2: " << v.flavor2()
            << " c2^dag: " << v.c_dagger_2() << " c2: " << v.c_2() << "\t" << v.abs_w();
}

std::ostream& operator<<(std::ostream& os, c_or_cdagger const& c) {
  return os << c.flavor() << " " << c.s() << " " << c.t();
}

std::ostream& operator<<(std::ostream& os, simple_hist const& h) {
  for (unsigned i = 0; i < h.size(); ++i) os << i << "\t" << h[i] << std::endl;
  return os;
}

std::ostream& operator<<(std::ostream& os, inverse_m_matrix const& m) {
  os << m.matrix() << std::endl << "creators: ";
  for (auto const& c : m.creators()) os << c << "\t";
  os << std::endl << "annihils: ";
  for (auto const& c : m.annihilators()) os << c << "\t";
  return os << std::endl;
}
