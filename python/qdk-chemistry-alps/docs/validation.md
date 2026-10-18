<!--
Copyright (c) 2026 Tobias Wolf
SPDX-License-Identifier: MIT

Adapted from skilledwolf/alpsx (private repository), branch
agent/v3-execution-and-api-audit at d37e51227,
python/qdk-chemistry-alps/docs/validation.md ("What the numbers establish"
and "Known limits"). Changes: the checks of this repository's suite; alpsx's
dated records (MPS bridge, Q# dynamics, historical benchmarks) are replaced
by a record of this repository's run.
-->

# Validation

The tests compare supported workflows against independent references: PySCF FCI, MACIS, `pyalps.ed`, direct small-system matrix calculations and QDK Chemistry's qubit mapper. They check energies, coefficients, observables and state and operator mappings; energy agreement alone cannot detect every convention error.

## Check your checkout

After [installing the solver stack](install.md), run from the **ALPS repository root**:

```bash
python -m pytest python/qdk-chemistry-alps --require-solver-stack -ra
```

`--require-solver-stack` makes a missing pyalps, qdk-chemistry or classic DMRG binding a failure instead of a skip. Inspect the skip report: the `RdmContainer` cases need a qdk-chemistry that provides it, and the `pyscf` cases need PySCF. Reinstall the SDK and pyalps after native changes before rerunning the suite.

## What the suite compares

- `alps_ed` against MACIS (energies and signed CI vectors) and PySCF FCI (H₂, H₄, and the sign-exact CI vector).
- `alps_classic_dmrg` against `alps_ed`, MACIS and PySCF FCI: energies, signed CI vectors, every 1-/2-RDM element and spin block, LiH in STO-3G with the default schedule; a truncated run's CI vector against its own energy; the measurements and an observable against `pyalps.ed`, and $\langle S^2\rangle$ against QDK's.
- The integral extraction, with the builder and with composite operators.
- The qubit export: the Fock-space spectrum for every encoding, and QDK's own mapper term by term.
- The state placement, through the energy, with a bit-reversed register as a negative control.
- The native Pauli route, against independently built encodings.
- Every script in [`examples/`](../examples/), each of which asserts its own physics (`tests/test_examples.py`): among them time-dependent DMRG and executed Q# circuits against exact propagation, and binary and Gray spin-1 encodings against an independent matrix reference.

## What the numbers establish

**Four-index Hamiltonians are checked beyond their energies.** Agreement of phase-aligned CI coefficients with PySCF tests fermion signs and Jordan–Wigner strings in the covered systems. A matching spectrum alone would not detect every basis convention error.

**The mode-ordering bridge is pinned.** ALPS uses interleaved ordering, qdk-chemistry and PySCF use blocked; `interleaved_to_blocked_parity` converts exactly, checked independently against both PySCF and MACIS. See [Conventions](conventions.md).

**DMRG observables are checked against FCI.** In the molecular cases, the RDMs of the final sweep match PySCF FCI element by element, spin block by spin block; these are the quantities an MCSCF driver consumes. The calculator's occupations and correlations are derived from these RDMs, and the suite compares them with `pyalps.ed` on the same Hamiltonian. Returning the RDMs as the wavefunction beyond the CI-vector limit needs a qdk-chemistry that provides `RdmContainer`.

**State mapping is checked through operator action.** A spectrum comparison is invariant under any permutation of the basis, so it cannot see an endianness or parity error. `ci_vector_to_qubit_statevector` places an ALPS state into the exported operator's register, and the operator's expectation in that state reproduces the ALPS energy within numerical tolerance. The test carries a negative control, a bit-reversed register that must fail the same check, on a model with a site-energy gradient, since a uniform chain's density profile is mirror-symmetric and would hide the bug.

**The chemistry-integral export route reuses QDK's fermion mappings.** For the supported untapered mappings, validated integral extraction feeds QDK Chemistry's Majorana engine and reproduces the Pauli terms of `QdkQubitMapper` on the same integrals. Tapering needs the separate QDK mapper workflow. The direct Pauli route uses the pyalps native compiler and an explicit layout; these are different contracts, as listed in [Capabilities](capabilities.md).

## Known limits

- **The calculator plugins target restricted electronic ground states.** They do not expose excited-state or finite-temperature solvers. The model-export workflow has a different contract; see [Capabilities and limits](capabilities.md).
- **RDM storage grows with K⁴.** Returning the full 2-RDM needs $K^4$ elements.
- **The calculator keeps no DMRG state.** Everything is measured in the final sweep; there is no continuation, and the calculator returns no matrix product state (`pyalps.classic_dmrg.solve_molecular` does, with `keep_mps`).
- **Performance records are workload-specific.** Dense synthetic integrals, sparse lattice interactions and molecular active spaces can have very different assembly and solver costs. Re-measure representative inputs on the checkout you plan to use.

## Record: 2026-10-18

With the commit that added the last of these pages, on macOS on Apple silicon, Python 3.13, PySCF 2.14:

| qdk-chemistry | qdk | NumPy | Result |
| --- | --- | --- | --- |
| 2.2.1 (PyPI) | 1.33.1 | 2.3.5 | 103 passed, 7 skipped (the `RdmContainer` cases) |
| 2.1.0 local build with `RdmContainer` | 1.32.3 | 2.5.3 | 108 passed, 2 skipped (the checks that only apply without `RdmContainer`) |

Each run takes about 90 s, about 80 s of it in the examples. The slowest example is the spin-1 quench, at about 30 s, almost all of it in the Q# simulator; the convergence scan takes about 7 s, and `large_active_space.py` about 6 s and 270 MB at 14 orbitals. Both qdk versions provide the Q# simulator, so the quench examples ran in both. The documentation's outputs and figures come from this state.
