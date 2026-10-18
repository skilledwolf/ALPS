<!--
Copyright (c) 2026 Tobias Wolf
SPDX-License-Identifier: MIT

Adapted from skilledwolf/alpsx (private repository), branch
agent/v3-execution-and-api-audit at d37e51227,
python/qdk-chemistry-alps/docs/examples/index.md. Changes: the three examples
that run on pyalps (the MPS, spin-1 and Q# quench examples are not carried
over); classic block DMRG in the solver background; plain commands; the
figures regenerated on this repository.
-->

# Examples

Explore a model, inspect its observables, and compare representations or computational methods. All three examples run with released qdk-chemistry 2.2.1. They study quantum models on ordinary computers.

| Example | Scientific question or technique |
| --- | --- |
| [Hubbard energies and correlations](hubbard-reference.md) | How do ED and DMRG agree on energies, occupations and spin correlations? Start here for the calculator workflow. |
| [Custom spin model → QDK](custom-spin.md) | How do multispin interactions and complex matrix elements survive operator conversion? |
| [Spinless fermion pairing](fermion-pairing.md) | How can a number-changing Hamiltonian cross the boundary without chemistry integrals? |

## Run one

After [setup](../install.md), run from the **ALPS repository root**:

```bash
python python/qdk-chemistry-alps/examples/hubbard_reference.py
```

Add `--output DIRECTORY` (needs `matplotlib`) to save a PNG plot and JSON results. The JSON records the model inputs, numerical controls, dependency versions and the git revision. The figures on these pages were saved that way on 2026-10-18; timings describe those runs, excluding imports and plotting, and are not performance guarantees.

## Shared notation and solver background

Sites are numbered from zero, and chains have open boundaries unless stated otherwise. We set $\hbar=1$: spin operators are dimensionless and energies are quoted in the stated hopping or exchange unit. For spin $1/2$, the Pauli operators satisfy $X=2S^x$, $Y=2S^y$ and $Z=2S^z$; a Pauli magnetization therefore ranges from $-1$ to $+1$.

**Exact diagonalization (ED)** solves $H\lvert\psi_n\rangle=E_n\lvert\psi_n\rangle$ in a finite basis, sometimes restricted by a conserved particle number or magnetization. "Exact" means within that basis and to numerical precision, not a solution without finite-size or basis-cutoff effects. It gives a reference for these small examples, but the Hilbert-space dimension grows exponentially with the number of sites.

**Density-matrix renormalization group (DMRG)** builds the chain from blocks and keeps, at each step, the $m$ most probable states of the reduced density matrix of a block, sweeping the boundary back and forth to optimize them. The kept states (`max_bond_dim`) bound the entanglement the state can represent across a cut; the state is a matrix product state of that bond dimension. Neither a successful run nor agreement between two observables establishes convergence: check the kept states, sweeps and truncation, as described in the [DMRG convergence guide](../how-to/converging-dmrg.md).

For an introduction to DMRG, see [White, *Density matrix formulation for quantum renormalization groups*](https://doi.org/10.1103/PhysRevLett.69.2863) and [Schollwöck, *The density-matrix renormalization group in the age of matrix product states*](https://arxiv.org/abs/1008.3477).

Each script checks its numerical results before reporting success. [Validation](../validation.md) identifies the independent references and environment-specific skips; [Conventions](../conventions.md) explains basis order and fermionic phases. A dirty-tree flag in saved JSON means the recorded revision alone does not reproduce that run.
