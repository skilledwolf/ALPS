<!--
Copyright (c) 2026 Tobias Wolf
SPDX-License-Identifier: MIT

Adapted from skilledwolf/alpsx (private repository), branch
agent/v3-execution-and-api-audit at d37e51227,
python/qdk-chemistry-alps/docs/how-to/molecules.md. Changes: alps_classic_dmrg
replaces alps_dmrg; plain commands; the historical ED size remark is replaced
by the classic DMRG route.
-->

# Run a molecule

The calculators accept real, restricted, spin-symmetric electronic integrals `(h1, eri, e_core)` in chemist notation. Up and down spins share the same spatial integrals, but their electron counts may differ. This example obtains them from PySCF; see [Capabilities and limits](../capabilities.md) for the exact input contract.

## From PySCF to a ground state

After [setup](../install.md), run from the **ALPS repository root** (needs `pyscf`):

```bash
python python/qdk-chemistry-alps/examples/molecule.py
```

The script, [`examples/molecule.py`](../../examples/molecule.py), uses LiH at a bond length of 1.6 Å in the STO-3G basis, with four electrons in six spatial orbitals, and solves it with `abinitio.build_abinitio_hamiltonian` and `solve_ground_state` (ED in the $(N, S_z)$ sector). Energies are in hartree.

```text
orbitals        : 6
Hartree-Fock    : -7.8618647698
ALPS ED         : -7.8823243789
PySCF FCI       : -7.8823243789
difference      : 7.11e-15
correlation     : -0.020460
```

This run agrees with PySCF full configuration interaction (FCI) to numerical precision within the chosen orbital basis. The [validation suite](../validation.md) checks several systems, including CI coefficients as well as energies.

## The two conventions that matter

**Use an orthonormal one-particle basis.** Here that is the Hartree–Fock molecular-orbital basis; raw atomic-orbital integrals generally use a nonorthogonal basis and cannot be passed directly. `ao2mo.kernel` transforms the two-electron integrals, and `ao2mo.restore(1, ...)` expands them to the full four-index array the builder expects.

**`eri[p, q, r, s]` uses chemist notation**, the $(pq|rs)$ ordering, eight-fold symmetric for real orbitals. Physicist notation uses a different index ordering; mixing them changes the Hamiltonian and may not raise an error. See [Conventions](../conventions.md#two-electron-integrals).

## Going through qdk-chemistry instead

Continuing with the integrals from the script, wrap them in a QDK `Hamiltonian` and use the registry:

```python
import qdk_chemistry_alps  # registers the calculators
from qdk_chemistry.algorithms import create
from qdk_chemistry.data import (
    CanonicalFourCenterHamiltonianContainer,
    Hamiltonian,
    ModelOrbitals,
)

container = CanonicalFourCenterHamiltonianContainer(
    h1, eri.reshape(-1), ModelOrbitals(n_orbitals), e_core, h1
)
hamiltonian = Hamiltonian(container)
calc = create("multi_configuration_calculator", "alps_ed")
energy, wavefunction = calc.run(hamiltonian, 2, 2)
```

Both ALPS calculators accept this `Hamiltonian`. Select `"alps_classic_dmrg"`, or a QDK calculator such as `"macis_cas"`, to compare methods on the same input: all three give `-7.8823243789` here. The final `h1` argument is the inactive Fock matrix; it equals the one-body integrals here because all orbitals are active. Frozen-core active spaces need their effective integrals and inactive Fock matrix.

## Active spaces

You can provide an active space instead of the full orbital set. Use the effective integrals for that space, its electron counts and the corresponding frozen-core energy. ED cost grows with the determinant count of the sector; for larger spaces use `alps_classic_dmrg`, assess its [convergence](converging-dmrg.md), and mind the separate [CI-vector limit](large-active-spaces.md).

## Next

- [Converge a DMRG run](converging-dmrg.md) once the active space outgrows ED
- [Measure observables](observables.md) to get RDMs for an MCSCF driver
