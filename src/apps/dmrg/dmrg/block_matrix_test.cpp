// SPDX-License-Identifier: MIT
#define WITH_LAPACK
#include <dmtk/dmtk.h>
#include <stdexcept>

template<class T>
void check_lookup(int count, int axis) {
  using namespace dmtk;
  QN::set_qn_mask(1 << axis);
  BMatrix<T> matrix;
  for (int i = 0; i < count; ++i) {
    QN qn;
    qn[axis] = 2*i - 4;
    SubSpace space(qn, 0, 0);
    matrix.push_back(SubMatrix<T>(qn, space, space));
  }
  auto const& constant = matrix;
  // Probe both ends, gaps between blocks, and an ignored quantum number.
  for (int value = -6; value <= 2*count; ++value) {
    QN query;
    query[axis] = value;
    query[1-axis] = 17;
    int index = (value + 4) / 2;
    if (value % 2 || value < -4 || index >= count) index = -1;
    auto expected = index < 0 ? nullptr : matrix[index];
    if (constant.get_block_index(query) != index ||
        constant.get_block(query) != expected || matrix.get_block(query) != expected)
      throw std::runtime_error("DMRG block lookup returned the wrong block");
  }
}

int main() {
  dmtk::QN::init();
  dmtk::QN::add_qn_index("first");
  dmtk::QN::add_qn_index("second");
  for (int count : {0, 1, 2, 3, 8, 17}) for (int axis : {0, 1}) {
    check_lookup<double>(count, axis);
    check_lookup<std::complex<double>>(count, axis);
  }
}
