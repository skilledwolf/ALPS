<!--
Copyright (c) 2026 Tobias Wolf
SPDX-License-Identifier: MIT

Adapted from skilledwolf/alpsx (private repository), branch
agent/v3-execution-and-api-audit at d37e51227,
python/qdk-chemistry-alps/docs/capabilities.md. Changes: the capabilities of
this repository (alps_classic_dmrg, pyalps.models export routes); the
prototype-scope notes, the MPS dynamics workflows and the dated
upstream review are replaced by this repository's state.
-->

# Capabilities and limits

Choose the route by **input representation and intended calculation**. The ED/DMRG calculator contract is narrower than the set of ALPS models that can be exported as qubit operators. The [export guide](how-to/hardware-export.md) provides code for the second route.

## QDK electronic Hamiltonian → ALPS calculation

| Capability | Contract |
| --- | --- |
| `alps_ed` / `alps_classic_dmrg` | Ground states of supported real, restricted electronic Hamiltonians, with explicit alpha/beta electron counts. "Spin-symmetric" means the same spatial integrals for both spins; the counts may differ. `alps_classic_dmrg` rejects integrals that break the real-orbital eightfold symmetry. |
| DMRG measurements | Local observables, correlations and reduced density matrices of the final sweep, derived from the RDMs; further operators through `observables` before the run. A full spatial 2-RDM has $K^4$ entries for $K$ orbitals. |
| Returned wavefunction | CI coefficients and determinants while the CI vector fits `ci_memory_limit_gb`. An RDM-only return needs a qdk-chemistry build providing `RdmContainer`; released qdk-chemistry 2.2.1 lacks it. |

These plugins do not accept arbitrary `QubitOperator`s or expose every native ALPS solver. The calculator keeps no DMRG state after the run: there is no continuation, and the calculator returns no matrix product state (`pyalps.classic_dmrg.solve_molecular(..., keep_mps=True)` does). Start with [Hubbard energies and correlations](examples/hubbard-reference.md).

## ALPS model → QDK qubit operator

| ALPS input | Export route and encoding | Main limit |
| --- | --- | --- |
| Spin-½ / hard-core bosons | `interop.hamiltonian_to_qubit_operator`: Pauli decomposition on two-level nonfermionic sites | Numeric coefficients; includes long-range and multispin products |
| Real, particle-conserving, spin-symmetric spinful fermions | `interop.fermionic_hamiltonian_to_qubit_operator`: restricted chemistry integrals → QDK Jordan–Wigner, Bravyi–Kitaev, BK-tree, parity or compatible untapered `MajoranaMapping` | Recognized ladder operators and real, eightfold-symmetric chemist integrals; only representable one-/two-body interactions; no tapering in this function |
| Strict spinless/spinful fermions, including complex hopping, pairing, spin mixing and higher-body terms | `pyalps.models.to_pauli_operator` → `interop.pauli_operator_to_qubit_operator`: native Jordan–Wigner | Dimension-2/4 binary occupation bases; Hermitian, even-parity Hamiltonians by default; local ALPS basis layout rather than chemistry blocked order |
| Positive integer/half-integer spins and finite bosonic spaces, including different local dimensions | The same native route, with binary or Gray encoding | Explicit finite bases and occupation cutoffs; binary is the default; bounded local matrix and Pauli expansion |

On the native route, export observables through `models.compile_observable` and the operator form of `to_pauli_operator`, reusing the Hamiltonian's encoding. Odd fermionic and non-Hermitian observables are supported on the corresponding spaces once `require_even` or `require_hermitian` is switched off. Keep the native layout alongside the QDK container: Pauli terms alone do not identify the physical basis or state map. All export routes have numerical tolerances; coefficient pruning can alter the represented operator.

Try the [custom spin export](examples/custom-spin.md), the [spin-1 encodings](examples/spin-one-anisotropy.md) or the [fermion pairing chain](examples/fermion-pairing.md).

## Other workflows

| Workflow | Current scope |
| --- | --- |
| ALPS lattice ↔ QDK `LatticeGraph` | Weighted adjacency conversion, with the geometry limitations listed below |
| Driven Hamiltonian | Supported static components can be exported separately; there is no adapter for the complete time-dependent protocol. |
| Time evolution | `pyalps.classic_dmrg.evolve` runs time-dependent DMRG on pyalps models; this package has no QDK interface for it. The [quench examples](examples/spin-quench.md) compare it with Q# product-formula circuits written from the exported Pauli terms; that circuit code is in the examples, not in this package. |

Constrained fermionic t–J sites and mixtures of higher-dimensional spins or bosons with fermions need additional encoding validation; general export rejects them. Automatic boson cutoffs, MPS-to-circuit preparation and a general protocol compiler are not supplied. Operator export by itself does not prepare a state or run a circuit, and RDMs do not provide coherent state amplitudes.

## Model catalogs versus expressiveness

ALPS describes models in two ways. The XML model and lattice catalogs drive the ALPS applications, and `pyalps.models.CatalogModel` compiles them for `pyalps.ed` and `pyalps.classic_dmrg`. `pyalps.models` also builds lattices and models programmatically, with ordered operator products similar in approach to QDK's programmatic interfaces: a product can involve four distinct sites or orbitals without being restricted to a graph edge, which is how `abinitio.build_abinitio_hamiltonian` represents general four-index two-body interactions. The supported QDK calculator inputs still have the restricted-integral contract listed above.

QDK Chemistry's named lattice-model constructors are Hückel, Hubbard, Pariser–Parr–Pople, anisotropic spin-½ Heisenberg and Ising. The Heisenberg constructor already covers cases that ALPS names separately, including XY and XXZ. These factories accept nonuniform parameters; arbitrary adjacency graphs extend beyond the built-in lattice shapes. See the [upstream model guide](https://microsoft.github.io/qdk-chemistry/user/comprehensive/model_hamiltonians.html).

The pyalps factories cover Heisenberg, Ising, XXZ, XY, transverse-field Ising, single-ion anisotropy, biquadratic exchange, J₁–J₂, Hubbard, interacting spinless fermions, t–J, Kondo lattice, hard-core bosons and Bose–Hubbard. The local operator algebra also supports higher spins, finite bosonic cutoffs, custom local matrices and mixed site types. The lattice layer keeps coordinates, unit cells, site and bond types and boundary information. A solver still has its own restrictions on which of these models it can execute.

The named catalogs are not the limits of either framework. QDK's Pauli representation can express any Hamiltonian on a finite qubit register. Higher-dimensional ALPS sites can be encoded into qubits, but that needs an explicit encoding and may expand the number of terms. Infinite bosonic spaces need a finite cutoff.

## Conventions and information that need care

- **Coupling normalization:** ALPS spin operators obey `S = σ/2` for spin ½. QDK model factories using Pauli matrices can therefore need different coupling and field coefficients. Export the operator instead of assuming identically named parameters match.
- **Basis and state mapping:** equal spectra do not establish correct qubit ordering or fermionic phases. Compare mapped states and observables as well. See [Conventions](conventions.md).
- **Geometry:** the lattice converter keeps adjacency, not coordinates, unit cells, site types or boundary-wrap metadata. Its assignment into an adjacency matrix does not accumulate parallel bonds, and it does not preserve or infer QDK term groupings. Export the expanded Hamiltonian when the complete interaction is what matters.
- **Solver metadata:** a Pauli sum alone does not transfer ALPS conserved-sector requests, measurements, convergence settings or a prepared state. A full-space qubit calculation and an ALPS sector calculation must be compared in matching physical sectors.
- **Scaling:** operator export avoids a full state vector. The chemistry route allocates a four-index integral array; the native Jordan–Wigner route avoids it, enforces a term-expansion bound and reports the discarded Pauli weight. Pauli expansion can still grow rapidly.

## Compatibility

Checked on 2026-10-18 against released qdk-chemistry 2.2.1 (qdk 1.33.1, NumPy 2.3.5) and a local qdk-chemistry 2.1.0 build providing `RdmContainer` (qdk 1.32.3, NumPy 2.5.3), both on Python 3.13 with PySCF 2.14. With 2.2.1, only the `RdmContainer` paths are unavailable. The [validation record](validation.md) lists the results.
