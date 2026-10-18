# Copyright (c) 2026 Tobias Wolf
# SPDX-License-Identifier: MIT
#
# Adapted from skilledwolf/alpsx (private repository), branch
# agent/v3-execution-and-api-audit at d37e51227,
# python/qdk-chemistry-alps/examples/large_active_space.py. Changes:
# alps_classic_dmrg replaces alps_dmrg: under the default
# wavefunction_container="auto" it returns an RdmContainer because the CI
# vector would exceed ci_memory_limit_gb; metadata key max_bond_dim.
"""Get a usable wavefunction past the CI-vector limit.

Requires a qdk-chemistry build with RdmContainer (released 2.2.1 lacks
it). At 14 orbitals half filled there are 11.8 million determinants:
about 94 MB of float64 amplitudes, but about 12 GB as qdk-chemistry
determinant objects, beyond the default ci_memory_limit_gb of 4. Under
the default wavefunction_container="auto" the calculator then returns the
active-space RDMs instead.
"""

import numpy as np
from qdk_chemistry.algorithms import create
from qdk_chemistry.data import LatticeGraph
from qdk_chemistry.utils.model_hamiltonians import create_hubbard_hamiltonian

import qdk_chemistry_alps  # noqa: F401  importing registers alps_ed and alps_classic_dmrg

n_orbitals = 14
hamiltonian = create_hubbard_hamiltonian(
    LatticeGraph.chain(n_orbitals), epsilon=0.0, t=1.0, U=4.0
)

calc = create("multi_configuration_calculator", "alps_classic_dmrg")
calc.settings().set("max_bond_dim", 64)
calc.settings().set("sweeps", 8)
energy, wavefunction = calc.run(hamiltonian, 7, 7)

print(f"container : {wavefunction.get_container_type()}")
print(f"energy    : {energy:.6f}")

container = wavefunction.get_container()
print(f"method    : {container.get_method()}")
print(f"kept states : {container.get_metadata()['max_bond_dim']}")

# The stored RDMs reproduce the energy through the integrals.
dm1 = np.asarray(container.get_active_one_rdm_spin_traced())
dm2 = np.asarray(container.get_active_two_rdm_spin_traced())
dm2 = dm2.reshape(n_orbitals, n_orbitals, n_orbitals, n_orbitals)

h1 = np.asarray(hamiltonian.get_one_body_integrals()[0])
two = hamiltonian.get_two_body_integrals()
eri = np.asarray(two[0] if isinstance(two, (tuple, list)) else two)
eri = eri.reshape(n_orbitals, n_orbitals, n_orbitals, n_orbitals)

from_rdms = (
    np.einsum("pq,pq->", h1, dm1)
    + 0.5 * np.einsum("pqrs,pqrs->", eri, dm2)
    + float(hamiltonian.get_core_energy())
)
print(f"energy from stored RDMs : {abs(from_rdms - energy):.2e}  (difference)")

# Natural occupations. Model bases carry no spin axis, so read the
# symmetry-blocked accessor rather than the spin-resolved pair.
occupations = wavefunction.active_orbital_occupations()
(label,) = occupations.extents()[0].keys()
natural = np.asarray(occupations.block([label])).ravel()
print(f"natural occupations sum : {natural.sum():.6f}  (= electron count)")

assert wavefunction.get_container_type() == "rdm"
assert abs(from_rdms - energy) < 1e-7
assert abs(natural.sum() - 14) < 1e-7
