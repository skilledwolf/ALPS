// SPDX-License-Identifier: MIT
#define WITH_LAPACK
#include <dmtk/dmtk.h>
#include <array>
#include <complex>
#include <stdexcept>

namespace {
using scalar = std::complex<double>;
using dense_vector = std::array<scalar, 16>;
struct action { int site; int shift; };

scalar element(int op, int column, int row, bool fermion) {
  if (fermion) return scalar(0.7 + 0.1 * op, 0.2);
  return scalar(0.2 * (1 + op + column + 2 * row), 0.1 * (op + column - row));
}

dense_vector apply(const dense_vector& input, action op, int index, bool fermion) {
  dense_vector result{};
  for (int state = 0; state < 16; ++state) {
    int bit = 1 << (4 - op.site), col = (state & bit) ? 1 : 0;
    for (int row = 0; row < 2; ++row) {
      if (fermion && row != col + op.shift) continue;
      int sign = 1;
      if (fermion)
        for (int site = 1; site < op.site; ++site)
          if (state & (1 << (4-site))) sign = -sign;
      result[(state & ~bit) | (row ? bit : 0)] +=
        double(sign) * element(index, col, row, fermion) * input[state];
    }
  }
  return result;
}

void check(std::vector<action> actions, bool fermion, bool hc, bool condensed = false) {
  using namespace dmtk;
  QN::init();
  QN::add_qn_index("N", true);
  QN::set_qn_mask(1);
  PackedBasis basis;
  if (fermion) {
    basis.push_back(SubSpace(QN(0), 0, 0));
    basis.push_back(SubSpace(QN(1), 1, 1));
  } else {
    basis.push_back(SubSpace(QN(0), 0, 1));
  }
  VectorState<scalar> input(basis, basis, basis, basis, QN(0), 0), output(input), cached(input);
  dense_vector original;
  for (int i = 0; i < 16; ++i) input[i] = original[i] = scalar(0.1 * (i+1), -0.03*i);
  output = cached = scalar(0);
  std::vector<BasicOp<scalar>> ops;
  for (unsigned i = 0; i < actions.size(); ++i) {
    BasicOp<scalar> op("test", 0);
    op.dqn = QN(fermion ? actions[i].shift : 0);
    op.set_fermion(fermion);
    op.resize(basis);
    for (auto& block : op) {
      for (unsigned col = 0; col < block.cols(); ++col)
        for (unsigned row = 0; row < block.rows(); ++row)
          block(col,row) = element(i, col, row, fermion);
    }
    ops.push_back(op);
  }
  scalar coefficient(0.7, -0.2);
  std::vector<ProductTerm<scalar>> terms;
  auto m = [&](int i) { return actions[i].site; };
  switch (ops.size()) {
  case 1:
    product(ops[0], input, output, m(0), coefficient, hc);
    terms = get_product_terms(ops[0], input, cached, m(0), coefficient, hc);
    break;
  case 2:
    product(ops[0], ops[1], input, output, m(0), m(1), coefficient, hc, static_cast<DMTKglobals<scalar>*>(nullptr), condensed);
    terms = get_product_terms(ops[0], ops[1], input, cached, m(0), m(1), coefficient, hc);
    break;
  case 3:
    product(ops[0], ops[1], ops[2], input, output, m(0), m(1), m(2), coefficient, hc);
    terms = get_product_terms(ops[0], ops[1], ops[2], input, cached, m(0), m(1), m(2), coefficient, hc);
    break;
  case 4:
    product(ops[0], ops[1], ops[2], ops[3], input, output, m(0), m(1), m(2), m(3), coefficient, hc);
    terms = get_product_terms(ops[0], ops[1], ops[2], ops[3], input, cached, m(0), m(1), m(2), m(3), coefficient, hc);
    break;
  }
  for (auto const& term : terms)
    product_term(term, input, cached, ops.size()==1 ? MASK_PRODUCT_DEFAULT : MASK_PRODUCT_DEFAULT|MASK_PRODUCT_HC);
  // Form a dense reference by applying local operators, including Jordan-Wigner
  // parity, to basis vectors. Adjoint handling is independent of the product code.
  dense_vector expected{};
  for (int col = 0; col < 16; ++col) {
    dense_vector column{}; column[col] = 1;
    for (int i = int(actions.size())-1; i >= 0; --i)
      column = apply(column, actions[i], i, fermion);
    for (int row = 0; row < 16; ++row) {
      scalar value = coefficient * column[row];
      if (!(actions.size()==1 && hc)) expected[row] += value * original[col];
      if (hc) expected[col] += std::conj(value) * original[row];
    }
  }
  for (int i = 0; i < 16; ++i) {
    if (std::abs(output[i] - expected[i]) > 1e-11 || std::abs(cached[i]-expected[i]) > 1e-11) {
      std::cerr << "component " << i << ": direct=" << output[i]
                << ", cached=" << cached[i] << ", expected=" << expected[i] << '\n';
      throw std::runtime_error("operator product differs from dense reference: arity=" +
        std::to_string(ops.size()) + " fermion=" + std::to_string(fermion) +
        " hc=" + std::to_string(hc));
    }
  }
}
}
int main() {
  for (bool hc : {false, true}) {
    for (int i=1;i<=4;++i) {
      check({{i,0}},false,hc);
      check({{i,1}},true,hc);
      for (int j=1;j<=4;++j) {
        check({{i,0},{j,0}},false,hc);
        // The condensed path requires a contiguous two-site slice.
        if (i==3 && j==4) check({{i,0},{j,0}},false,hc,true);
        if (i<=j) check({{i,1},{j,-1}},true,hc);
      }
    }
    for (auto sites : {std::vector<int>{1,2,3}, {1,2,4}, {1,3,4}, {2,3,4}, {4,2,1}, {1,2,3,4}, {4,3,2,1}}) {
      std::vector<action> actions;
      for (int site : sites) actions.push_back({site,0});
      check(actions,false,hc);
    }
    check({{1,1},{2,-1},{3,1},{4,-1}},true,hc);
  }
}
