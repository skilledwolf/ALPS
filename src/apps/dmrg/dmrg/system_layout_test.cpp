// SPDX-License-Identifier: MIT
#define WITH_LAPACK
#include <dmtk/dmtk.h>
#include <array>
#include <stdexcept>

namespace {
using namespace dmtk;

struct test_system : System<double> {
  void warmup(bool value) { _in_warmup = value; }
};

void check(std::array<int, 4> lengths, int total, std::vector<int> normal,
           std::vector<int> inward, std::vector<int> outward) {
  for (int flags = 0; flags < 8; ++flags) for (int singles = 0; singles < 16; ++singles) {
    test_system system;
    system.h.lattice() = Lattice(total, OBC);
    system.warmup(flags & 1);
    system.set_grow_symmetric(flags & 2);
    system.set_grow_outward(flags & 4);
    std::array<Block<double>, 4> blocks;
    system._b1 = &blocks[0]; system._b2 = &blocks[1];
    system._b3 = &blocks[2]; system._b4 = &blocks[3];
    bool growing = (flags & 1) && (flags & 2);
    auto const& owners = growing ? ((flags & 4) ? outward : inward) : normal;
    std::vector<BasicOp<double> const*> expected(total, nullptr);
    for (int b = 0; b < 4; ++b) {
      blocks[b].set_lattice(lengths[b] ? Lattice(lengths[b], OBC) : Lattice());
      blocks[b].set_single_site(singles & (1 << b));
      int local = 0;
      for (int site = 0; site < total; ++site) if (owners[site] == b + 1) {
        blocks[b].push_back(BasicOp<double>("probe", blocks[b].single_site() ? local : site));
        expected[site] = &blocks[b].back();
        ++local;
      }
    }
    for (int site = -2; site < total + 2; ++site) {
      auto owner = site >= 0 && site < total ? owners[site] : BLOCK_NONE;
      auto op = site >= 0 && site < total ? expected[site] : nullptr;
      if (system.block(site) != owner || system(BasicOp<double>("probe", site)) != op ||
          system(BasicOp<double>("missing", site)) != nullptr)
        throw std::runtime_error("DMRG block/site mapping changed");
    }
  }
}
}

int main() {
  QN::init();
  check({2,1,1,3}, 7, {1,1,2,3,4,4,4}, {1,1,2,3,4,4,4}, {1,1,2,3,4,4,4});
  check({2,1,1,3}, 12,
        {1,1,2,3,4,4,4,0,0,0,0,0}, {1,1,2,0,0,0,0,0,3,4,4,4},
        {0,0,0,1,1,2,3,4,4,4,0,0});
  check({2,1,1,3}, 13,
        {1,1,2,3,4,4,4,0,0,0,0,0,0}, {1,1,2,0,0,0,0,0,0,3,4,4,4},
        {0,0,0,1,1,2,3,4,4,4,0,0,0});
  check({0,1,1,0}, 4, {2,3,0,0}, {2,0,0,3}, {0,2,3,0});
}
