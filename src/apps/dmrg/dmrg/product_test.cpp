// SPDX-License-Identifier: MIT
#define WITH_LAPACK
#include <dmtk/dmtk.h>
#include <array>
#include <complex>
#include <stdexcept>

namespace {
using scalar = std::complex<double>;
using dense_vector = std::vector<scalar>;
struct action { int site; int shift; };

scalar element(int op, int column, int row, bool fermion) {
  if (fermion) return scalar(0.7 + 0.1 * op, 0.2);
  return scalar(0.2 * (1 + op + column + 2 * row), 0.1 * (op + column - row));
}

dense_vector apply(const dense_vector& input, action op, int index, bool fermion,
                   const std::array<int, 4>& dimensions) {
  dense_vector result(input.size());
  int stride = 1;
  for (int site = op.site; site < 4; ++site) stride *= dimensions[site];
  for (unsigned state = 0; state < input.size(); ++state) {
    int col = (state / stride) % dimensions[op.site-1];
    for (int row = 0; row < dimensions[op.site-1]; ++row) {
      if (fermion && row != col + op.shift) continue;
      int sign = 1;
      if (fermion)
        for (int site = 1; site < op.site; ++site)
          if (state & (1 << (4-site))) sign = -sign;
      result[state + (row-col)*stride] +=
        double(sign) * element(index, col, row, fermion) * input[state];
    }
  }
  return result;
}

void check(std::vector<action> actions, bool fermion, bool hc, bool condensed = false,
           std::array<int, 4> dimensions = {2,2,2,2}) {
  using namespace dmtk;
  QN::init();
  QN::add_qn_index("N", true);
  QN::set_qn_mask(1);
  std::array<PackedBasis, 4> bases;
  int size = 1;
  for (int site = 0; site < 4; ++site) {
    size *= dimensions[site];
    if (fermion) {
      bases[site].push_back(SubSpace(QN(0), 0, 0));
      bases[site].push_back(SubSpace(QN(1), 1, 1));
    } else {
      bases[site].push_back(SubSpace(QN(0), 0, dimensions[site]-1));
    }
  }
  VectorState<scalar> input(bases[0], bases[1], bases[2], bases[3], QN(0), 0), output(input), cached(input);
  dense_vector original(size);
  for (int i = 0; i < size; ++i) input[i] = original[i] = scalar(0.1 * (i+1), -0.03*i);
  output = cached = scalar(0);
  std::vector<BasicOp<scalar>> ops;
  for (unsigned i = 0; i < actions.size(); ++i) {
    BasicOp<scalar> op("test", 0);
    op.dqn = QN(fermion ? actions[i].shift : 0);
    op.set_fermion(fermion);
    op.resize(bases[actions[i].site-1]);
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
  dense_vector expected(size), normal(size), adjoint(size);
  for (int col = 0; col < size; ++col) {
    dense_vector column(size); column[col] = 1;
    for (int i = int(actions.size())-1; i >= 0; --i)
      column = apply(column, actions[i], i, fermion, dimensions);
    for (int row = 0; row < size; ++row) {
      scalar value = coefficient * column[row];
      normal[row] += value * original[col];
      adjoint[col] += std::conj(value) * original[row];
      if (!(actions.size()==1 && hc)) expected[row] += value * original[col];
      if (hc) expected[col] += std::conj(value) * original[row];
    }
  }
  for (int i = 0; i < size; ++i) {
    if (std::abs(output[i] - expected[i]) > 1e-11 || std::abs(cached[i]-expected[i]) > 1e-11) {
      std::cerr << "component " << i << ": direct=" << output[i]
                << ", cached=" << cached[i] << ", expected=" << expected[i] << '\n';
      throw std::runtime_error("operator product differs from dense reference: arity=" +
        std::to_string(ops.size()) + " fermion=" + std::to_string(fermion) +
        " hc=" + std::to_string(hc));
    }
  }
  if (ops.size() == 2) {
    DMTKglobals<scalar> scratch;
    for (auto* globals : {static_cast<DMTKglobals<scalar>*>(nullptr), &scratch})
      for (int contribution : {0, int(MASK_PRODUCT_DEFAULT), int(MASK_PRODUCT_HC),
                               MASK_PRODUCT_DEFAULT | MASK_PRODUCT_HC}) {
        const scalar initial(0.4, -0.3);
        cached = initial;
        for (const auto& term : terms)
          product_term2(term, input, cached, contribution, globals, condensed);
        for (int i = 0; i < size; ++i) {
          scalar selected = initial;
          if (contribution & MASK_PRODUCT_DEFAULT) selected += normal[i];
          if (hc && (contribution & MASK_PRODUCT_HC)) selected += adjoint[i];
          if (std::abs(cached[i] - selected) > 1e-11)
            throw std::runtime_error("two-operator contribution differs from dense reference");
        }
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
        check({{i,0},{j,0}},false,hc,false,{2,3,4,2});
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
