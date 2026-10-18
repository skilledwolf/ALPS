<!--
Copyright (c) 2026 Tobias Wolf
SPDX-License-Identifier: MIT

Adapted from skilledwolf/alpsx (private repository), branch
agent/v3-execution-and-api-audit at d37e51227,
python/qdk-chemistry-alps/docs/examples/index.md. Changes: classic block DMRG
and its time-dependent variant in the solver background, in place of MPS
DMRG and evolution; plain commands; the figures regenerated on this
repository; the presentation plan and concept-diagram links are not carried
over.
-->

# Examples

Explore a model, inspect its observables, and compare representations or computational methods. All seven examples run with released qdk-chemistry 2.2.1; the two quenches also need the `qdk` package's Q# simulator. ALPS reference calculations and Q# circuit simulations both study quantum models on ordinary computers.

| Example | Scientific question or technique |
| --- | --- |
| [Hubbard energies and correlations](hubbard-reference.md) | How do ED and DMRG agree on energies, occupations and spin correlations? Start here for the calculator workflow. |
| [Custom spin model → QDK](custom-spin.md) | How do multispin interactions and complex matrix elements survive operator conversion? |
| [Spin-chain quench](spin-quench.md) | How do exact propagation, time-dependent DMRG and simulated circuit measurements compare? |
| [Spin-1 anisotropy](spin-one-anisotropy.md) | How does single-ion anisotropy redistribute local spin populations, and how are three levels encoded? |
| [Spin-1 dynamics and leakage](spin-one-quench.md) | How do binary and Gray circuit approximations trade accuracy, leakage and gate count? Read the anisotropy example first. |
| [Spinless fermion pairing](fermion-pairing.md) | How can a number-changing Hamiltonian cross the boundary without chemistry integrals? |
| [Fourteen-site RDM observables](rdm-observables.md) | What can be measured from a DMRG state without reconstructing all CI amplitudes? This uses ALPS directly. |

## Run one

After [setup](../install.md), run from the **ALPS repository root**:

```bash
python python/qdk-chemistry-alps/examples/hubbard_reference.py
```

Add `--output DIRECTORY` (needs `matplotlib`) to save a PNG plot and JSON results. The JSON records the model inputs, numerical controls, dependency versions and the git revision; the quench JSON also holds the Q# program and the sampled counts. The figures on these pages were saved that way on 2026-10-18; timings describe those runs, excluding imports and plotting, and are not performance guarantees.

## Shared notation and solver background

Sites are numbered from zero, and chains have open boundaries unless stated otherwise. We set $\hbar=1$: spin operators are dimensionless, energies are quoted in the stated hopping or exchange unit, and time is measured in the corresponding inverse energy unit. The dynamics pages use $\tau$ for time to distinguish it from the Hubbard hopping $t$. For spin $1/2$, the Pauli operators satisfy $X=2S^x$, $Y=2S^y$ and $Z=2S^z$; a Pauli magnetization therefore ranges from $-1$ to $+1$.

**Exact diagonalization (ED)** solves $H\lvert\psi_n\rangle=E_n\lvert\psi_n\rangle$ in a finite basis, sometimes restricted by a conserved particle number or magnetization. "Exact" means within that basis and to numerical precision, not a solution without finite-size or basis-cutoff effects. It gives a reference for these small examples, but the Hilbert-space dimension grows exponentially with the number of sites.

**Density-matrix renormalization group (DMRG)** builds the chain from blocks and keeps, at each step, the $m$ most probable states of the reduced density matrix of a block, sweeping the boundary back and forth to optimize them. The kept states (`max_bond_dim`) bound the entanglement the state can represent across a cut; the state is a matrix product state of that bond dimension. Neither a successful run nor agreement between two observables establishes convergence: check the kept states, sweeps and truncation, as described in the [DMRG convergence guide](../how-to/converging-dmrg.md).

**Time-dependent DMRG (tDMRG)**, `pyalps.classic_dmrg.evolve`, evolves a state in the same block basis: second-order Trotter gates on neighbouring sites are applied while sweeping, and the kept basis follows the state ([White and Feiguin](https://doi.org/10.1103/PhysRevLett.93.076401)). It adds a time-step error to the truncation error; the spin-chain quench varies both.

**QDK export and Q# evolution** are different steps. Export expresses a Hamiltonian as a sum of Pauli operators on qubits. The quench examples then approximate $e^{-\mathrm{i}H\tau}$ by a sequence of Pauli rotations and execute that circuit in QDK's software simulator. A smaller circuit time step reduces the product-formula error; more measurement shots reduce sampling uncertainty. Neither change fixes too few kept DMRG states. Each example explains which of these quantities its figure actually compares.

For an introduction to DMRG, see [White, *Density matrix formulation for quantum renormalization groups*](https://doi.org/10.1103/PhysRevLett.69.2863) and [Schollwöck, *The density-matrix renormalization group in the age of matrix product states*](https://arxiv.org/abs/1008.3477).

Each script checks its numerical results before reporting success. [Validation](../validation.md) identifies the independent references and environment-specific skips; [Conventions](../conventions.md) explains basis order and fermionic phases. A dirty-tree flag in saved JSON means the recorded revision alone does not reproduce that run.
