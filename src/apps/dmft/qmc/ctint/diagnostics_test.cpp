// SPDX-License-Identifier: MIT
#include "types.hpp"
#include <iostream>
#include <sstream>
#include <stdexcept>

int main() {
  creator c(1, 2, 0.25, 0);
  vertex v(0, 1, 2, 3, 1, 2, 4, 5, 0.5);
  vertex_array vertices{v};
  simple_hist histogram(2);
  histogram[1] = 7;
  inverse_m_matrix matrix;
  matrix.matrix().resize(1, 1);
  matrix.matrix()(0, 0) = 2;
  matrix.creators().emplace_back(1, 2, 0.25, 0);
  matrix.annihilators().emplace_back(0, 3, 0.5, 0);
  std::ostringstream out, console, expected;
  auto* previous = std::cout.rdbuf(console.rdbuf());
  bool chained = &(out << c) == &out;
  out << '\n' << v << '\n' << vertices << histogram << matrix;
  std::cout.rdbuf(previous);
  std::string vertex_text = "z1: 0 c1^dag: 2 c1: 3 z2: 1 c2^dag: 4 c2: 5\t0.5\n";
  expected << "1 2 0.25\n" << vertex_text << vertex_text << "0\t0\n1\t7\n"
           << matrix.matrix() << std::endl
           << "creators: 1 2 0.25\t\nannihils: 0 3 0.5\t\n";
  if (!chained || !console.str().empty() || out.str() != expected.str())
    throw std::runtime_error("CTINT diagnostics ignored their destination stream");
}
