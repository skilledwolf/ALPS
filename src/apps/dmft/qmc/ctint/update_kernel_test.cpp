// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
#include "update_kernel.hpp"
#include <alps/numeric/matrix/algorithms.hpp>
#include <map>
#include <stdexcept>
#include <string>

namespace {
void require(bool condition) {
  if (!condition) throw std::runtime_error("CT-INT kernel contract failed");
}

struct matrix_state {
  alps::numeric::matrix<double> values;
  std::vector<int> creation, annihilation;
  std::vector<double> shifts;
  auto& matrix() { return values; }
  auto& creators() { return creation; }
  auto& annihilators() { return annihilation; }
  auto& alpha() { return shifts; }
};

struct matrix_run {
  std::vector<matrix_state> M{1};
  double green0_spline(int row, int col) const {
    const double green[3][3] = {{2, .5, .25}, {.5, 3, .4}, {.25, .4, 4}};
    return green[row][col];
  }
  void check_inverse(unsigned size) {
    alps::numeric::matrix<double> direct(size, size);
    for (unsigned i = 0; i < size; ++i)
      for (unsigned j = 0; j < size; ++j)
        direct(i, j) = green0_spline(M[0].creation[i], M[0].annihilation[j]) + (i == j ? .1 : 0);
    auto expected = inverse(direct);
    for (unsigned i = 0; i < size; ++i)
      for (unsigned j = 0; j < size; ++j)
        require(std::abs(expected(i, j) - M[0].values(i, j)) < 1.e-12);
  }
};

struct samples : std::vector<double> {
  void operator<<(double value) { push_back(value); }
};

struct workflow_run {
  std::vector<int> vertices;
  unsigned max_order = 3;
  double sign = 1, weight = 0;
  std::map<std::string, samples> measurements;
  const std::vector<double> random_values{.1, .2, .9, .0, .1, .1, .95, .9, .1};
  unsigned draws = 0, rejected = 0;
  double ctint_uniform() { return random_values.at(draws++); }
  double try_add() { return -.5; }
  double try_remove(unsigned vertex) { require(vertex == 0); return -.5; }
  void perform_add() { vertices.push_back(0); }
  void perform_remove(unsigned) { vertices.pop_back(); }
  void reject_add() { ++rejected; }
  void reject_remove() { ++rejected; }
};
}

int main() {
  using kernel = alps::ctint::update_kernel;
  matrix_run matrices;
  for (unsigned i = 0; i < 3; ++i) {
    matrices.M[0].creation.push_back(i);
    matrices.M[0].annihilation.push_back(i);
    matrices.M[0].shifts.push_back(.1);
    double proposal = kernel::fastupdate_up(matrices, 0, true);
    require(num_rows(matrices.M[0].values) == i);
    require(proposal == kernel::fastupdate_up(matrices, 0, false));
    matrices.check_inverse(i + 1);
  }
  double proposal = kernel::fastupdate_down(matrices, 1, 0, true);
  require(num_rows(matrices.M[0].values) == 3);
  require(proposal == kernel::fastupdate_down(matrices, 1, 0, false));
  require(num_rows(matrices.M[0].values) == 2);
  matrices.check_inverse(2);

  workflow_run run;
  kernel::interaction_expansion_step(run); // accepted insertion
  require(run.vertices.size() == 1 && run.sign == -1 && run.draws == 2);
  kernel::interaction_expansion_step(run); // accepted removal
  require(run.vertices.empty() && run.sign == 1 && run.draws == 5);
  kernel::interaction_expansion_step(run); // rejected insertion
  require(run.vertices.empty() && run.rejected == 1 && run.draws == 7);
  kernel::interaction_expansion_step(run); // removal from an empty configuration
  require(run.draws == 8 && run.weight == -.5);
  run.max_order = 0;
  kernel::interaction_expansion_step(run); // insertion at the order limit
  require(run.draws == 9 && run.vertices.empty());
  require(run.measurements["VertexInsertion"] == std::vector<double>({1, 0}));
  require(run.measurements["VertexRemoval"] == std::vector<double>({1}));
}
