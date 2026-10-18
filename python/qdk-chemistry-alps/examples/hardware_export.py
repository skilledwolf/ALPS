# Copyright (c) 2026 Tobias Wolf
# SPDX-License-Identifier: MIT
#
# Adapted from skilledwolf/alpsx (private repository), branch
# agent/v3-execution-and-api-audit at d37e51227,
# python/qdk-chemistry-alps/examples/hardware_export.py. Changes: pyalps
# names (pyalps.models, models.lattice.chain).
"""Export one fermionic Hamiltonian as a QDK qubit operator.

The ALPS and QDK representations must describe the same quantum model.
This example checks their full spectra by numerical diagonalization on an
ordinary computer; it does not construct a circuit or run on a device.
Spectral agreement alone does not establish a correct state/basis map.
"""

import numpy as np
from pyalps import models

from qdk_chemistry_alps import interop

# A fermionic model built on the ALPS side.
H = models.hubbard(models.lattice.chain(3, bc="open"), t=1.0, U=4.0, mu=0.7)

for encoding in ["jordan-wigner", "bravyi-kitaev", "parity"]:
    qubit_op = interop.fermionic_hamiltonian_to_qubit_operator(H, encoding)

    # Compare the full spectra, not just the ground state: the operators
    # must agree on every state, constant term included.
    fock = np.sort(np.linalg.eigvalsh(H.to_sparse_matrix().toarray()))
    as_qubits = interop.qubit_operator_to_hamiltonian(qubit_op)
    qubit = np.sort(np.linalg.eigvalsh(as_qubits.to_sparse_matrix().toarray()))

    terms = len(qubit_op.to_dict()["pauli_strings"])
    print(
        f"{encoding:18s} {qubit_op.num_qubits} qubits, {terms:2d} Pauli terms, "
        f"spectrum agrees to {np.max(np.abs(fock - qubit)):.1e}"
    )
    assert np.max(np.abs(fock - qubit)) < 1e-12

# Metadata travels with the operator, so downstream code knows the layout.
meta = interop.fermionic_hamiltonian_to_qubit_operator(H).to_dict()
print(f"encoding={meta['encoding']!r} mode_order={meta['fermion_mode_order']!r}")
