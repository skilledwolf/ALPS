# Copyright (c) 2026 Tobias Wolf
# SPDX-License-Identifier: MIT
#
# Adapted from skilledwolf/alpsx (private repository), branch
# agent/v3-execution-and-api-audit at d37e51227,
# python/qdk-chemistry-alps/examples/quickstart.py. Changes: alps_classic_dmrg
# (classic block DMRG) replaces alps_dmrg (MPS DMRG).
"""Solve a Hubbard chain with both backends and compare them.

The shortest end-to-end path: build a model with qdk-chemistry's own
factory, solve it two ways, and check they agree.
"""

from qdk_chemistry.algorithms import create
from qdk_chemistry.data import LatticeGraph
from qdk_chemistry.utils.model_hamiltonians import create_hubbard_hamiltonian

import qdk_chemistry_alps  # noqa: F401  importing registers alps_ed and alps_classic_dmrg

# A six-site Hubbard chain at half filling.
hamiltonian = create_hubbard_hamiltonian(
    LatticeGraph.chain(6), epsilon=0.0, t=1.0, U=4.0
)

# Exact diagonalization in the full fixed-particle-number sector.
ed = create("multi_configuration_calculator", "alps_ed")
energy_ed, wavefunction = ed.run(hamiltonian, 3, 3)

# Classic block DMRG, checked against ED here.
dmrg = create("multi_configuration_calculator", "alps_classic_dmrg")
dmrg.settings().set("max_bond_dim", 128)
energy_dmrg, _ = dmrg.run(hamiltonian, 3, 3)

print(f"exact diagonalization : {energy_ed:.10f}")
print(f"DMRG                  : {energy_dmrg:.10f}")
print(f"difference            : {abs(energy_ed - energy_dmrg):.2e}")
print(f"determinants returned : {len(wavefunction.get_coefficients())}")

assert abs(energy_ed - energy_dmrg) < 1e-8
