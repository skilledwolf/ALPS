<!--
Copyright (c) 2026 Tobias Wolf
SPDX-License-Identifier: MIT

Adapted from skilledwolf/alpsx (private repository), branch
agent/v3-execution-and-api-audit at d37e51227,
python/qdk-chemistry-alps/docs/how-to/hardware-export.md. Changes: pyalps names
(models.to_pauli_operator, compile_observable for observables,
models.lattice); the comparison figure is not carried over; output
re-measured.
-->

# Export a qubit operator

Build a supported model in ALPS and export its Hamiltonian as a QDK `QubitOperator`. The result describes a Hamiltonian on a qubit register. It does not prepare a state, choose an algorithm or execute a circuit. The [quench examples](../examples/spin-quench.md) add explicit circuit construction and **software simulation**; running on quantum hardware would require further preparation, compilation and execution steps.

The native route supports finite spin and boson spaces, with binary and Gray encoding, and strict spinless and spinful fermions, including pairing. The separate chemistry route below supports several fermion encodings for the real, spin-symmetric, particle-conserving class representable as restricted chemistry integrals. See [Capabilities and limits](../capabilities.md) for the precise boundaries.

## Exporting a higher-spin model

```python
from pyalps import models
from qdk_chemistry_alps import interop

H = models.heisenberg(models.lattice.chain(4, bc="open"), S=1)
encoding = models.make_qubit_encoding(H, scheme="binary")
encoded = models.to_pauli_operator(H, encoding=encoding)
qubit_op = interop.pauli_operator_to_qubit_operator(encoded)
sz = models.to_pauli_operator(
    models.compile_observable(H, H.site_for(0).Sz(0)), models.physical_basis(H), encoding=encoding
)
initial_bits = encoding.encode_basis_bits([0, 1, 2, 1])
```

Keep `encoding` alongside the QDK container and use it for observables and state preparation. Four spin-1 sites occupy eight qubits with some unused local states; see [higher spins](../conventions.md#higher-spins-and-finite-local-spaces). The [anisotropy example](../examples/spin-one-anisotropy.md) explains the physics; the [quench example](../examples/spin-one-quench.md) executes Q# circuits and checks leakage. The same native route handles explicit finite boson cutoffs and general supported fermions; see the [pairing example](../examples/fermion-pairing.md).

## Fermions through chemistry integrals

After [setup](../install.md), run from the **ALPS repository root**:

```bash
python python/qdk-chemistry-alps/examples/hardware_export.py
```

[`examples/hardware_export.py`](../../examples/hardware_export.py) builds a three-site Hubbard chain in ALPS, exports it with three encodings and compares the full spectra:

```python
H = models.hubbard(models.lattice.chain(3, bc="open"), t=1.0, U=4.0, mu=0.7)
qubit_op = interop.fermionic_hamiltonian_to_qubit_operator(H, "bravyi-kitaev")
```

```text
jordan-wigner      6 qubits, 18 Pauli terms, spectrum agrees to 8.0e-15
bravyi-kitaev      6 qubits, 18 Pauli terms, spectrum agrees to 8.0e-15
parity             6 qubits, 18 Pauli terms, spectrum agrees to 9.8e-15
encoding='jordan-wigner' mode_order='blocked'
```

Available encodings: `"jordan-wigner"` (default), `"bravyi-kitaev"`, `"bravyi-kitaev-tree"`, `"parity"`, or any untapered `MajoranaMapping`. They are unitarily equivalent (same spectrum, different Pauli weight and connectivity), so the choice is about what your device and algorithm prefer, not about the physics.

For tapered encodings such as symmetry-conserving Bravyi–Kitaev, map the qdk-chemistry Hamiltonian with `QdkQubitMapper`, which applies tapering as a post-pass. Passing a tapered mapping here raises rather than silently dropping the tapering.

## What the spectrum check establishes

The full spectrum checks every energy level and the constant offset, which is stronger than checking the ground-state energy alone. However, any basis permutation preserves the spectrum: this check cannot establish correct qubit order or fermionic state phases. The [validation suite](../validation.md#what-the-numbers-establish) also maps states, compares expectation values and uses an incorrect-ordering negative control. New encodings need both kinds of check.

## Spin-½ convenience adapter

For spin-½ models, the convenience adapter remains available:

```python
from pyalps import models
from qdk_chemistry_alps import interop

H_spin = models.heisenberg(models.lattice.chain(4, bc="open"), S=0.5)
qubit_op = interop.hamiltonian_to_qubit_operator(H_spin)
```

This takes ALPS Hamiltonians on spin-½ sites and decomposes each term in the Pauli basis. Fermionic sites are rejected: they need an encoding, which the native route or the chemistry-integral route above provides.

## Closing the comparison loop

For a meaningful comparison, keep the physical Hamiltonian, state and observables consistent across representations. Building the input once avoids accidentally changing couplings between branches. Independent constructions are still valuable when validating the exporter itself. This snippet prepares two representations; it does not run a solver or a circuit:

```python
from pyalps import models
from qdk_chemistry_alps import interop

# 1. Build once.
H = models.hubbard(models.lattice.chain(6, bc="open"), t=1.0, U=4.0)

# 2. Export a QDK qubit operator.
qubit_op = interop.fermionic_hamiltonian_to_qubit_operator(H, "jordan-wigner")

# 3. Extract integrals from the same H for an electronic solver.
h1, eri, e_core = interop.hamiltonian_to_integrals(H)
```

For the blocked Jordan–Wigner route, `ci_vector_to_qubit_statevector` maps ALPS CI amplitudes into the exported register: the operator's expectation in that state reproduces the ALPS energy, and the tests check exactly that, with a bit-reversed register as a negative control. It does not perform the additional state transformation needed for Bravyi–Kitaev or parity encodings, and it is not a circuit state-preparation method.

The [Hubbard example](../examples/hubbard-reference.md) shows the solver call and [observable measurements](observables.md). For a separately executed dynamics comparison, continue with the [Q# quench](../examples/spin-quench.md), which matches initial states and observation times as well as the Hamiltonian. Correlations test more of the state than energy alone.

## Conventions to carry across the boundary

**Blocked fermion mode order on the chemistry route.** Modes 0…K−1 are alpha and K…2K−1 are beta. Only Jordan–Wigner directly assigns each occupation to its matching qubit; Bravyi–Kitaev and parity use different encoded bits. The `fermion_mode_order="blocked"` metadata describes the input mode layout, not a universal qubit-occupation rule.

**Little-endian Pauli strings.** The rightmost character acts on qubit 0.

**Explicit native layout.** The higher-spin and direct Jordan–Wigner route keeps the ordered ALPS local states with site 0 least significant. Use its `encoding` and state helpers; chemistry blocked or interleaved metadata does not describe it.

**Compare matrices after applying the state map.** The spin-½ convenience route needs bit reversal because ALPS and QDK use opposite site significance. Higher-spin and fermionic encodings need their own maps, including phases where applicable. Spectral agreement alone does not establish a correct mapping.

See [Conventions](../conventions.md) for the full set.

## Next

- [Measure observables](observables.md) to build the ED or DMRG reference
- [Converge a DMRG run](converging-dmrg.md) so the reference deserves the name
