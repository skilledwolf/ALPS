// SPDX-License-Identifier: MIT
#include <alps/model.h>
#include <alps/model/blochbasisstates.h>
#include <functional>
#include <numeric>
#include <stdexcept>

using constraints_type = std::vector<std::pair<std::string, alps::half_integer<short>>>;

void check(alps::basis_states_descriptor<short> const& basis, constraints_type const& constraints) {
  std::vector<std::vector<short>> expected;
  std::vector<short> state(basis.size());
  // Independent exhaustive enumeration, without partial-bound pruning.
  std::function<void(size_t)> enumerate = [&](size_t site) {
    if (site < basis.size()) {
      for (short index=0; index<basis[site].size(); ++index) {
        state[site] = index;
        enumerate(site+1);
      }
      return;
    }
    for (auto const& constraint : constraints) {
      alps::half_integer<short> total;
      for (size_t i=0; i<basis.size(); ++i)
        total += get_quantumnumber(basis[i][state[i]], constraint.first, basis.get_site_basis(i));
      if (total != constraint.second) return;
    }
    expected.push_back(state);
  };
  if (!basis.empty()) enumerate(0);
  alps::basis_states<short> actual(basis, constraints);
  if (actual.size() != expected.size() || !std::equal(actual.begin(), actual.end(), expected.begin()))
    throw std::runtime_error("constrained basis differs from exhaustive enumeration");

  // An identity-only translation leaves every accepted state as its own orbit.
  std::vector<size_t> identity(basis.size());
  std::iota(identity.begin(), identity.end(), 0);
  alps::bloch_basis_states<short>::translation_type translations{{1., identity}};
  alps::bloch_basis_states<short> bloch(basis, translations, constraints);
  if (bloch.full_list() != expected || bloch.size() != expected.size())
    throw std::runtime_error("Bloch consumer changed basis enumeration order");
  for (size_t i=0; i<expected.size(); ++i) {
    auto entry = bloch.index_and_phase(expected[i]);
    if (entry.first != i || entry.second != std::complex<double>(1) ||
        std::abs(bloch.normalization(i) - 1/std::sqrt(double(basis.size()))) > 1e-14)
      throw std::runtime_error("identity Bloch projection changed");
  }
}

int main(int argc, char** argv) {
  if (argc != 2) return 1;
  check(alps::basis_states_descriptor<short>(), {});
  for (std::string model : {"spin", "spinless fermions"}) {
    alps::Parameters p;
    p["MODEL_LIBRARY"] = std::string(argv[1]) + "/models.xml";
    p["LATTICE_LIBRARY"] = std::string(argv[1]) + "/lattices.xml";
    p["LATTICE"] = "open chain lattice";
    p["L"] = 4;
    p["MODEL"] = model;
    alps::ModelLibrary models(p);
    alps::graph_helper<> lattice(p);
    auto ham = models.get_hamiltonian(model);
    p.copy_undefined(ham.default_parameters());
    ham.set_parameters(p);
    alps::basis_states_descriptor<short> basis(ham.basis(), lattice.graph());
    check(basis, {});
    std::string quantum_number = model == "spin" ? "Sz" : "N";
    for (int total : {0,1,2,999})
      check(basis, {{quantum_number, alps::half_integer<short>(total)}});
    check(basis, {{quantum_number, alps::half_integer<short>(0)},
                  {quantum_number, alps::half_integer<short>(1)}});
  }
}
