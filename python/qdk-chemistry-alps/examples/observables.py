# Copyright (c) 2026 Tobias Wolf
# SPDX-License-Identifier: MIT
#
# Adapted from skilledwolf/alpsx (private repository), branch
# agent/v3-execution-and-api-audit at d37e51227,
# python/qdk-chemistry-alps/examples/observables.py. Changes:
# alps_classic_dmrg replaces alps_dmrg; the measurements come from the run's
# RDMs, not from an MPS contraction.
"""Measure observables on a converged DMRG state.

The classic engine measures every 1-/2-RDM element in its final sweep;
occupations and correlations are derived from them, so they stay
affordable beyond the CI limit. Solver convergence must be checked
separately.
"""

import numpy as np
from qdk_chemistry.algorithms import create
from qdk_chemistry.data import LatticeGraph
from qdk_chemistry.utils.model_hamiltonians import create_hubbard_hamiltonian

import qdk_chemistry_alps  # noqa: F401  importing registers alps_ed and alps_classic_dmrg

n_orbitals = 6
hamiltonian = create_hubbard_hamiltonian(
    LatticeGraph.chain(n_orbitals), epsilon=0.0, t=1.0, U=4.0
)

calc = create("multi_configuration_calculator", "alps_classic_dmrg")
calc.settings().set("max_bond_dim", 128)
energy, _ = calc.run(hamiltonian, 3, 3)

# Orbital occupations, and the density-density and spin-spin correlations.
occupations = np.asarray(calc.measure_local("n")).real
density = np.asarray(calc.measure_correlations("n")).real
spin = np.asarray(calc.spin_spin_correlations()).real

print(f"energy            : {energy:.6f}")
print(f"occupations       : {np.round(occupations, 4)}")
print(f"trace             : {occupations.sum():.6f}  (= electron count)")

neighbour = [spin[i, i + 1] for i in range(n_orbitals - 1)]
print(f"<S_i.S_j> nearest : {np.round(neighbour, 4)}")
print(f"mean              : {np.mean(neighbour):.4f}  (negative = antiferromagnetic)")

# measure_local and measure_correlations are read from these RDMs, so this
# checks the index conventions rather than the state.
dm1, dm2 = calc.make_rdms()
dm1 = np.asarray(dm1).real
dm2 = np.asarray(dm2).real
index = np.arange(n_orbitals)
density_from_rdms = dm2[index, index][:, index, index] + np.diag(np.diag(dm1))

print(f"occupations vs 1-RDM diagonal : {np.max(np.abs(occupations - np.diag(dm1))):.2e}")
print(f"<n_i n_j> vs 2-RDM            : {np.max(np.abs(density - density_from_rdms)):.2e}")

assert abs(occupations.sum() - 6) < 1e-8
assert np.mean(neighbour) < -0.1
assert np.max(np.abs(occupations - np.diag(dm1))) < 1e-9
assert np.max(np.abs(density - density_from_rdms)) < 1e-9
