<!--
Copyright (c) 2026 Tobias Wolf
SPDX-License-Identifier: MIT

Adapted from skilledwolf/alpsx (private repository), branch
agent/v3-execution-and-api-audit at d37e51227,
python/qdk-chemistry-alps/docs/reference/api.md. Changes: a written summary of
this package's entry points replaces the mkdocstrings pages; the docstrings
hold the full parameter descriptions.
-->

# API reference

Start with the [quickstart](../quickstart.md) for a complete calculation; use [Capabilities](../capabilities.md) to choose a supported input and [Conventions](../conventions.md) to interpret states and observables. The docstrings hold the full parameter descriptions (`help(qdk_chemistry_alps.interop)`).

| Task | Entry point |
| --- | --- |
| Solve an electronic Hamiltonian | Import `qdk_chemistry_alps`, then select `alps_ed` or `alps_classic_dmrg` with QDK's `create` registry |
| Measure the last DMRG run | `make_rdms`, `measure_local`, `measure_correlations`, `spin_spin_correlations`, `last_result` |
| Convert a model or encoding | [`interop`](#qdk_chemistry_alpsinterop); native encodings belong to `pyalps.models` |
| Build directly from electronic integrals | [`abinitio`](#qdk_chemistry_alpsabinitio) |

The calculator call is `energy, wavefunction = calc.run(hamiltonian, n_alpha, n_beta)`. Its input is a QDK electronic `Hamiltonian`; it is not a pyalps Hamiltonian or a QDK `QubitOperator`. Settings are described [separately](settings.md).

## `qdk_chemistry_alps`

Importing the package registers the calculators. It exports `AlpsEdCalculator`, `AlpsClassicDmrgCalculator`, `build_abinitio_hamiltonian`, the `interop` and `abinitio` modules, `load` and `__version__`. `AlpsClassicDmrgCalculator` is registered only when pyalps has its `classic_dmrg` bindings.

## Calculators

**`AlpsEdCalculator`** (`"alps_ed"`): exact diagonalization of the integrals as a `pyalps.models` Hamiltonian, one `SpinfulFermion()` site per orbital, in the $(N, S_z)$ sector with `pyalps.ed`. Returns the CI vector as a `StateVectorContainer`.

**`AlpsClassicDmrgCalculator`** (`"alps_classic_dmrg"`): classic block DMRG of restricted real integrals through `pyalps.classic_dmrg.solve_molecular`, which builds the Hamiltonian, sweeps and measures every 1-/2-RDM element by spin block in the final sweep and, for small active spaces, contracts the CI vector from the engine's stored rotations. Returns a `StateVectorContainer` with the CI vector and the run's RDMs, or an `RdmContainer`. On the calculator:

| Member | Returns |
| --- | --- |
| `make_rdms()` | Spin-traced `(dm1, dm2)` in the PySCF `make_rdm12` convention |
| `measure_local(op)` | $\langle O_i\rangle$ per orbital for `n`, `n_up`, `n_dn`, `Sz`, `n_updn` |
| `measure_correlations(op_i, op_j=None)` | $C_{ij}=\langle O_i^\dagger O_j\rangle$ for number-type operators, or one of `c_up`, `c_dn`, `Sp`, `Sm` |
| `spin_spin_correlations()` | $\langle\mathbf S_i\cdot\mathbf S_j\rangle$; its sum is $\langle S^2\rangle$ |
| `observables` | Set before `run()`: a dict of `pyalps.models` term sums measured in the final sweep |
| `last_result` | The `pyalps.classic_dmrg` result of the last run |

They read the most recent `run()` and raise `RuntimeError` before one.

## `qdk_chemistry_alps.interop`

| Function | Converts |
| --- | --- |
| `lattice_graph_to_alps(graph)` | QDK `LatticeGraph` → pyalps custom `Lattice` and its bond-type weights |
| `alps_lattice_to_lattice_graph(lattice, type_weights=None)` | pyalps `Lattice` → QDK `LatticeGraph` (adjacency only) |
| `qubit_operator_to_hamiltonian(qubit_op, lattice=None, conserved_qns=())` | QDK `QubitOperator` → pyalps `Hamiltonian` on spin-½ sites |
| `hamiltonian_to_qubit_operator(H, tol=1e-12)` | pyalps `Hamiltonian` on two-level bosonic sites → `QubitOperator` |
| `pauli_operator_to_qubit_operator(operator)` | native `pyalps.models.PauliOperator` → `QubitOperator` |
| `hamiltonian_to_integrals(H, tol=1e-10)` | particle-conserving, spin-symmetric strict-fermion `Hamiltonian` → `(h1, eri, e_core)` |
| `fermionic_hamiltonian_to_qubit_operator(H, encoding="jordan-wigner", ...)` | the same class of `Hamiltonian` → `QubitOperator` through QDK's fermion mappings |
| `ci_vector_to_qubit_statevector(basis_states, coefficients, n_orbitals)` | ALPS CI vector → amplitudes on the blocked Jordan–Wigner register |
| `rcm_ordering(graph)` | Reverse Cuthill–McKee site ordering of a `LatticeGraph` |

## `qdk_chemistry_alps.abinitio`

| Function | Does |
| --- | --- |
| `build_abinitio_hamiltonian(h1, eri, *, tol=1e-12)` | A pyalps `Hamiltonian` of chemist-notation integrals, one `SpinfulFermion()` site per orbital, `N` and `Sz` conserved |
| `solve_ground_state(H, n_alpha, n_beta, *, k=1)` | `pyalps.ed` eigenstates in the $(N, S_z)$ sector: `eigenvalues`, `eigenvectors`, `problem` |
| `alps_state_to_occupations(state_row)` | One ALPS basis row → `(alpha_occ, beta_occ)` |
| `interleaved_to_blocked_parity(alpha_occ, beta_occ)` | The sign between the ALPS and PySCF mode orders |

## `qdk_chemistry_alps.plugin`

`load()` registers the ALPS algorithms with qdk-chemistry; it is idempotent. qdk-chemistry builds with plugin discovery also load the `qdk_chemistry.plugins` entry point `alps`.
