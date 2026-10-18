# Copyright (c) 2026 Tobias Wolf
# SPDX-License-Identifier: MIT
#
# Adapted from skilledwolf/alpsx (private repository), branch
# agent/v3-execution-and-api-audit at d37e51227,
# python/qdk-chemistry-alps/examples/interop.py. Changes: pyalps names
# (pyalps.models, models.lattice.chain, ed.ground_state for
# alps.ed.eigenstates).
"""Move models between qdk-chemistry and ALPS representations.

Both packages can express the same lattice systems in different terms, so
these adapters translate at the boundary instead of duplicating either
side's model factories.
"""

import numpy as np
from pyalps import ed, models
from qdk_chemistry.data import LatticeGraph
from qdk_chemistry.utils.model_hamiltonians import create_heisenberg_hamiltonian

from qdk_chemistry_alps import interop

# --- a spin model built in qdk-chemistry, solved in ALPS -------------------
qubit_op = create_heisenberg_hamiltonian(
    LatticeGraph.chain(6), jx=1.0, jy=1.0, jz=1.0, hz=0.3
)
H = interop.qubit_operator_to_hamiltonian(qubit_op)
energy = ed.ground_state(H).energy
print(f"Heisenberg chain ground state : {energy:.6f}")

# Check that this supported spin operator round-trips within tolerance.
returned = interop.hamiltonian_to_qubit_operator(H)
original_terms = dict(zip(
    qubit_op.to_dict()["pauli_strings"],
    qubit_op.to_dict()["coefficients"]["real"],
    strict=True,
))
returned_terms = dict(zip(
    returned.to_dict()["pauli_strings"],
    returned.to_dict()["coefficients"]["real"],
    strict=True,
))
worst = max(abs(original_terms[k] - returned_terms[k]) for k in original_terms)
print(f"round trip worst coefficient  : {worst:.1e}")

# --- a fermionic model built in ALPS, handed to qdk-chemistry solvers ------
H_fermi = models.hubbard(models.lattice.chain(4, bc="open"), t=1.0, U=4.0)
h1, eri, e_core = interop.hamiltonian_to_integrals(H_fermi)
print(f"extracted integrals           : h1 {h1.shape}, eri {eri.shape}")

# --- orbital ordering for a two-dimensional lattice -----------------------
# DMRG maps a lattice onto a chain, so its cost depends on how far apart
# connected sites end up. A lattice arriving in arbitrary order is the case
# worth fixing -- the natural row-major order of a square lattice is already
# near-optimal, so shuffle first to show what RCM actually buys.
shuffled = np.random.default_rng(0).permutation(12)
adjacency = np.asarray(LatticeGraph.square(3, 4).adjacency_matrix(), dtype=float)
adjacency = adjacency[shuffled][:, shuffled]
order = interop.rcm_ordering(LatticeGraph.from_dense_matrix(adjacency))


def bandwidth(matrix):
    rows, cols = np.nonzero(matrix)
    return int(np.max(np.abs(rows - cols))) if rows.size else 0


print(
    f"3x4 lattice bandwidth         : {bandwidth(adjacency)} -> "
    f"{bandwidth(adjacency[order][:, order])} after reordering"
)

assert worst < 1e-12
assert h1.shape == (4, 4)
assert bandwidth(adjacency[order][:, order]) < bandwidth(adjacency)
