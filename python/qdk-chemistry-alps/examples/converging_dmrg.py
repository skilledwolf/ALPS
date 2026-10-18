# Copyright (c) 2026 Tobias Wolf
# SPDX-License-Identifier: MIT
#
# Adapted from skilledwolf/alpsx (private repository), branch
# agent/v3-execution-and-api-audit at d37e51227,
# python/qdk-chemistry-alps/examples/converging_dmrg.py. Changes:
# alps_classic_dmrg replaces alps_dmrg; the RdmContainer setting is dropped,
# as the calculator returns this small problem's CI vector with released
# qdk-chemistry too.
"""Check that a DMRG run is converged, rather than assuming it.

This scan checks energy stability under an increasing number of kept
states. Also check sweeps, Hamiltonian tolerances and the observables of
interest.
"""

from qdk_chemistry.algorithms import create
from qdk_chemistry.data import LatticeGraph
from qdk_chemistry.utils.model_hamiltonians import create_hubbard_hamiltonian

import qdk_chemistry_alps  # noqa: F401  importing registers alps_ed and alps_classic_dmrg

hamiltonian = create_hubbard_hamiltonian(
    LatticeGraph.chain(8), epsilon=0.0, t=1.0, U=8.0
)

energies = {}
for bond_dimension in [4, 8, 16, 32, 64, 128]:
    calc = create("multi_configuration_calculator", "alps_classic_dmrg")
    calc.settings().set("max_bond_dim", bond_dimension)
    energies[bond_dimension], _ = calc.run(hamiltonian, 4, 4)

best = min(energies.values())
print(f"{'bond dim':>9}  {'energy':>14}  {'above best':>12}")
for bond_dimension, energy in energies.items():
    print(f"{bond_dimension:>9}  {energy:>14.9f}  {energy - best:>12.2e}")

# For this case, check that energy decreases as the variational space grows.
# Finite-sweep optimization need not reach the minimum in each space.
values = list(energies.values())
assert all(later <= earlier + 1e-9 for earlier, later in zip(values, values[1:], strict=False))

# Agreement checks this bond-dimension increment, not the exact-state error.
converged = abs(values[-1] - values[-2])
print(f"\nlast two differ by {converged:.2e} -- converged at this precision")
assert converged < 1e-6
