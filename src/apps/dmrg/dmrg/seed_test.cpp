// SPDX-License-Identifier: MIT
#define WITH_LAPACK
#include <dmtk/dmtk.h>
#include <stdexcept>
#include <type_traits>

namespace {
using namespace dmtk;
int case_number = 0;

template<class T>
struct seed_system : System<T> {
  seed_system(int length, int iteration, int direction, bool basic) {
    this->_lattice = Lattice(length, OBC);
    this->h = Hami<T>(this->_lattice);
    for (int i = 0; i < length; ++i) this->h.sites[i] = nullptr;
    this->set_name(("seed_test_" + std::to_string(++case_number)).c_str());
    this->iter = iteration;
    this->dir = direction;
    this->_use_basic_seed = basic;
    this->_in_warmup = false;
    this->set_qn_mask(0);
    this->_target.resize(0);
  }
  using System<T>::init_iteration;
};

template<class T>
void check(int direction, int iteration, int rho_index, int other_index,
           int gs_index, int position, int length = 8, bool basic = false, bool use_seed = true) {
  seed_system<T> system(length, iteration, direction, basic);
  Basis basis(1);
  basis.reorder();
  Basis pair(basis, basis);
  pair.reorder();
  Block<T> block(basis);
  block.set_lattice(Lattice(1, OBC));
  VectorState<T> state;
  state.set_qn_mask(QN(), 0);
  state.resize(basis, basis, basis, basis);
  state = T(3);
  system.gs = state;
  system.seed = state;
  system.seed = T(13);
  system._propagate_state.resize(1);
  system._propagate_state[0] = state;
  system._propagate_state[0] = T(5);
  T seed_expected = !use_seed ? T(13) : basic ? T(1) : T(3);
  T propagated_expected = T(5);
  if (use_seed && !basic && length != 4) {
    T r1 = T(2), r2 = T(7);
    if constexpr (std::is_same_v<T, std::complex<double>>) { r1 = T(2, 1); r2 = T(7, -2); }
    SubSpace space(QN(), 0, 0);
    SubMatrix<T> value(QN(), space, space);
    BMatrix<T> rho;
    value = r1;
    rho.push_back(value);
    system.write_rho(rho, pair, rho_index, position);
    *rho[0] = r2;
    system.write_rho(rho, pair, other_index, position == LEFT ? RIGHT : LEFT);
    system.write_gs(state, gs_index, position);
    // Poison the in-memory states to ensure init_iteration reads the checkpoint.
    system.gs = T(17);
    system._propagate_state[0] = T(19);
    seed_expected = conj(r1)*r2*T(3);
    propagated_expected = conj(r1)*r2*T(5);
  }
  system.init_iteration(block, block, block, block, use_seed, false);
  if (system.seed.size() != 1 || std::abs(system.seed[0] - seed_expected) > 1e-12 ||
      system._propagate_state.size() != 1 || system._propagate_state[0].size() != 1 ||
      std::abs(system._propagate_state[0][0] - propagated_expected) > 1e-12)
    throw std::runtime_error("DMRG seed or propagated state changed");
}

template<class T>
void check_seeds() {
  check<T>(RIGHT2LEFT, 1, 5, 2, 4, LEFT);
  check<T>(LEFT2RIGHT, 1, 5, 2, 4, RIGHT);
  check<T>(LEFT2RIGHT, 2, 2, 5, 1, LEFT);
  check<T>(RIGHT2LEFT, 2, 2, 5, 1, RIGHT);
  for (int direction : {LEFT2RIGHT, RIGHT2LEFT}) {
    check<T>(direction, 1, 0, 0, 0, LEFT, 4);
    check<T>(direction, 1, 0, 0, 0, LEFT, 4, true);
    check<T>(direction, 2, 0, 0, 0, LEFT, 8, true);
    check<T>(direction, 2, 0, 0, 0, LEFT, 8, false, false);
  }
}
}

int main() {
  QN::init();
  auto directory = boost::filesystem::temp_directory_path() / boost::filesystem::unique_path("alps-seed-%%%%-%%%%");
  boost::filesystem::create_directory(directory);
  tmp_files.set_temp_dir(directory.string().c_str());
  check_seeds<double>();
  check_seeds<std::complex<double>>();
  boost::filesystem::remove_all(directory);
}
