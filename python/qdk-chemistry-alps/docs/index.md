<!--
Copyright (c) 2026 Tobias Wolf
SPDX-License-Identifier: MIT

Adapted from skilledwolf/alpsx (private repository), branch
agent/v3-execution-and-api-audit at d37e51227,
python/qdk-chemistry-alps/docs/index.md and docs/how-to/index.md. Changes: the
calculators are alps_ed and alps_classic_dmrg on pyalps; the prototype-scope
notes, the MPS and Q# workflows and the development pages are not carried over.
-->

# qdk-chemistry-alps

Use ALPS solvers from **QDK Chemistry**, or export an ALPS model as a **QDK qubit operator**. This optional Python package lives in the ALPS repository at `python/qdk-chemistry-alps/` and runs on [pyalps](../../pyalps/README.md).

## Start here

1. [Install](install.md) the ALPS SDK, pyalps and this package.
2. [Run the quickstart](quickstart.md): solve a six-site Hubbard chain with exact diagonalization (ED) and classic block DMRG.
3. [Explore the examples](examples/index.md): Hubbard correlations against an independent reference, a custom spin model and a fermion pairing chain exported to QDK.

## Choose a workflow

**Solve a QDK electronic Hamiltonian.** Select `alps_ed` or `alps_classic_dmrg` through QDK's calculator registry. Both return an energy and a wavefunction; the DMRG calculator also keeps the run's reduced density matrices (RDMs) for [occupations and correlation functions](how-to/observables.md). Start with [Hubbard energies and correlations](examples/hubbard-reference.md) or [molecular integrals](how-to/molecules.md).

**Export an ALPS model to QDK.** Translate supported spin, finite-boson and fermionic operators into Pauli terms while keeping their basis and state maps. Circuit construction and execution are separate steps. Start with [qubit-operator export](how-to/hardware-export.md).

The calculator plugins accept restricted electronic ground-state problems; operator export supports a broader model class. [Capabilities and limits](capabilities.md) identifies the route for a given input. Returning an RDM-only QDK wavefunction needs a qdk-chemistry build that provides `RdmContainer`; everything else works with released qdk-chemistry 2.2.1.

## How-to guides

| Guide | Answers |
| --- | --- |
| [Solve a lattice model](how-to/lattice-models.md) | How do I get a ground state for a Hubbard or Heisenberg system? |
| [Run a molecule](how-to/molecules.md) | How do I go from PySCF integrals to an ALPS calculation? |
| [Measure observables](how-to/observables.md) | How do I get RDMs, occupations and correlation functions? |
| [Work past the CI limit](how-to/large-active-spaces.md) | What do I do when the active space is too big for a CI vector? |
| [Export a qubit operator](how-to/hardware-export.md) | How do I transfer a Hamiltonian and its conventions to QDK? |
| [Converge a DMRG run](how-to/converging-dmrg.md) | How do I know the number I got is right? |

## Reference

| Need | Read |
| --- | --- |
| A setting, function or convention | [Settings](reference/settings.md), [API](reference/api.md), [Conventions](conventions.md) |
| What runs where | [Capabilities and limits](capabilities.md) |
| Evidence behind the results | [Validation](validation.md) |
| How the figures were made | [Figures](figures/README.md) |

Every guide's code is a script in [`examples/`](../examples/), and `tests/test_examples.py` runs each of them, so the commands and outputs on these pages are checked with the test suite.
