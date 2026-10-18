<!--
Copyright (c) 2026 Tobias Wolf
SPDX-License-Identifier: MIT

Adapted from skilledwolf/alpsx (private repository), branch
agent/v3-execution-and-api-audit at d37e51227,
python/qdk-chemistry-alps/docs/quickstart.md. Changes: alps_classic_dmrg
replaces alps_dmrg; plain commands instead of the Pixi environment; output
re-measured on this repository.
-->

# Quickstart

Solve a six-site Hubbard chain with both ALPS backends through QDK Chemistry. This checks the installed solver workflow and introduces the inputs and results used throughout the guides.

## Run it

Complete [Installation](install.md) first. From the **ALPS repository root**:

```bash
python python/qdk-chemistry-alps/examples/quickstart.py
```

The script, [`examples/quickstart.py`](../examples/quickstart.py), builds the model with qdk-chemistry's own factory and compares the two ground-state energies:

```python
hamiltonian = create_hubbard_hamiltonian(LatticeGraph.chain(6), epsilon=0.0, t=1.0, U=4.0)

ed = create("multi_configuration_calculator", "alps_ed")
energy_ed, wavefunction = ed.run(hamiltonian, 3, 3)

dmrg = create("multi_configuration_calculator", "alps_classic_dmrg")
dmrg.settings().set("max_bond_dim", 128)
energy_dmrg, _ = dmrg.run(hamiltonian, 3, 3)
```

```text
exact diagonalization : -3.0925653195
DMRG                  : -3.0925653195
difference            : 1.64e-14
determinants returned : 400
```

Last digits may vary; the script requires agreement within `1e-8`.

## Understand the result

The chain has six spatial orbitals and three electrons of each spin, passed as `run(hamiltonian, 3, 3)`. The fixed-particle-number basis contains $\binom{6}{3}^2=400$ determinants. Energies are in units of the hopping, $t=1$; the onsite repulsion is $U=4$.

- **Registration:** `import qdk_chemistry_alps` makes `alps_ed` and `alps_classic_dmrg` available in QDK's calculator registry.
- **Return values:** each call returns `(energy, wavefunction)`. Here `wavefunction` is the ED result, with CI coefficients and determinants. The DMRG calculator returns its own CI vector for a problem this small and keeps the run's RDMs for later measurements.
- **Accuracy:** ED solves within the selected finite basis to numerical precision. DMRG optimizes a truncated block basis. Agreement checks this small case; larger or different models need their own [convergence checks](how-to/converging-dmrg.md).

Both plugins accept supported **restricted, spin-symmetric electronic Hamiltonians**. A `QubitOperator` is a different input type and cannot be passed to these calculators. See [Capabilities](capabilities.md) before changing the model class.

## Next

Follow [Hubbard energies and correlations](examples/hubbard-reference.md) for the physics and a plot, or [measure observables](how-to/observables.md) on the DMRG state. For orbital integrals from PySCF, continue with [Run a molecule](how-to/molecules.md).
