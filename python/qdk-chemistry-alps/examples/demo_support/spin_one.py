# Copyright (c) 2026 Tobias Wolf
# SPDX-License-Identifier: MIT
#
# Adapted from skilledwolf/alpsx (private repository), branch
# agent/v3-execution-and-api-audit at d37e51227,
# python/qdk-chemistry-alps/examples/demo_support/spin_one.py. Changes: pyalps
# names (models.Hamiltonian, models.Spin, models.lattice.chain,
# conserved_quantum_numbers); qubit_matrix is dropped for the package's
# pauli_matrix, which reads the same QDK container; encoding_layout stands in
# for alpsx's QubitEncoding.describe().
"""Four-site spin-1 models and small independent matrix references."""

import numpy as np
from pyalps import models

COUPLINGS = np.array([1., .8, 1.2])
FIELDS = np.array([.07, -.11, .16, -.03])


def embed(matrix, site, n=4):
    result = np.ones((1, 1), dtype=complex)
    for i in range(n):
        result = np.kron(result, matrix if i == site else np.eye(3))
    return result


def model(anisotropy, *, biquadratic=.18, transverse=0., y_field=0.):
    spin = models.Spin(S=1.)
    h = models.Hamiltonian(models.lattice.chain(4, bc="open"), spin,
                           conserved_quantum_numbers=["Sz"] if transverse == y_field == 0 else [])
    # Explicit angular momentum matrices, independent of ALPS operators.
    x = np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]]) / np.sqrt(2)
    y = np.array([[0, -1j, 0], [1j, 0, -1j], [0, 1j, 0]]) / np.sqrt(2)
    z = np.diag([1., 0., -1.])
    axes = [[embed(s, i) for i in range(4)] for s in (x, y, z)]
    reference = .37 * np.eye(81, dtype=complex)
    h.add(.37 * spin.I(0))
    for i, coupling in enumerate(COUPLINGS):
        expression = (spin.Sz(i)*spin.Sz(i+1) + .5*spin.Sp(i)*spin.Sm(i+1)
                      + .5*spin.Sm(i)*spin.Sp(i+1))
        h.add(coupling*expression + biquadratic*expression*expression)
        exchange = sum(axis[i] @ axis[i+1] for axis in axes)
        reference += coupling*exchange + biquadratic*exchange@exchange
    for i, field in enumerate(FIELDS):
        h.add(anisotropy*spin.Sz(i)*spin.Sz(i) + field*spin.Sz(i))
        reference += anisotropy*axes[2][i]@axes[2][i] + field*axes[2][i]
    if transverse:
        h.add(transverse*spin.Sx(1))
        reference += transverse*axes[0][1]
    if y_field:
        h.add(y_field*spin.Sy(2))
        reference += y_field*axes[1][2]
    return h, reference, axes


def reference_indices(scheme):
    """Independent declared mapping: physical site 0 is most significant."""
    indices = []
    for local in np.ndindex(3, 3, 3, 3):
        codes = [j if scheme == "binary" else j ^ (j >> 1) for j in local]
        indices.append(sum(code << (2*i) for i, code in enumerate(codes)))
    return np.asarray(indices)


def encoding_layout(encoding):
    """The register layout of a pyalps QubitEncoding, for the saved data."""
    return {"scheme": encoding.scheme.name, "num_qubits": encoding.num_qubits,
            "site_dimensions": list(encoding.site_dimensions), "site_qubits": list(encoding.site_qubits),
            "site_offsets": list(encoding.site_offsets)}
