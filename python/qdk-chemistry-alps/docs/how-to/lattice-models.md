<!--
Copyright (c) 2026 Tobias Wolf
SPDX-License-Identifier: MIT

Adapted from skilledwolf/alpsx (private repository), branch
agent/v3-execution-and-api-audit at d37e51227,
python/qdk-chemistry-alps/docs/how-to/lattice-models.md. Changes: pyalps names
(pyalps.models factories, models.lattice, pyalps.ed.ground_state); classic
block DMRG instead of an MPS; plain commands.
-->

# Solve a lattice model

Start from QDK Chemistry for its electronic-Hamiltonian solver workflow, or from ALPS for its broader set of local spaces and interactions. A spin model represented as a QDK `QubitOperator` takes the conversion route below; it is not an input to the electronic calculator plugins.

## From a qdk-chemistry factory

Build a six-site Hubbard chain, then solve its half-filled sector with ALPS:

```python
from qdk_chemistry.data import LatticeGraph
from qdk_chemistry.utils.model_hamiltonians import create_hubbard_hamiltonian

graph = LatticeGraph.chain(6)
hamiltonian = create_hubbard_hamiltonian(graph, epsilon=0.0, t=1.0, U=4.0)
```

`LatticeGraph` also provides `square`, `honeycomb`, `kagome` and `triangular` geometries. `LatticeGraph.from_dense_matrix` accepts a custom adjacency matrix.

```python
import qdk_chemistry_alps  # registers the calculators
from qdk_chemistry.algorithms import create

calc = create("multi_configuration_calculator", "alps_ed")
energy, wavefunction = calc.run(hamiltonian, 3, 3)
print(energy)
```

The two `3` arguments count spin-up and spin-down electrons. If you change the graph size or filling, adjust them accordingly. The [quickstart](../quickstart.md) adds a DMRG comparison, and the [Hubbard example](../examples/hubbard-reference.md) explains the physics and measures correlations.

## From an ALPS model factory

pyalps builds Hamiltonians programmatically: named factories (`heisenberg`, `xxz`, `hubbard`, `t_j`, `bose_hubbard`, ...) on `models.lattice` lattices, local spaces such as `Spin(S)`, `Boson(n_max)` or `SpinfulFermion()`, and explicit operator products. The XML model catalog is available too, through `models.CatalogModel`. For a spin chain:

```python
from pyalps import ed, models

spin = models.heisenberg(models.lattice.chain(6, bc="open"), J=1.0)
result = ed.ground_state(spin, sector={"Sz": 0})
print(result.energy)   # -2.4935771338879276
```

This uses spin operators, $H=J\sum_i\mathbf S_i\cdot\mathbf S_{i+1}$, with $S^\alpha=\sigma^\alpha/2$ for spin $1/2$. QDK's Heisenberg factory in the interop example below uses Pauli coefficients; equal numerical couplings therefore do not imply equal Hamiltonians. See [Conventions](../conventions.md). The [pyalps README](../../../pyalps/README.md#models-exact-diagonalization-and-classic-dmrg) covers the model layer, `pyalps.ed` and `pyalps.classic_dmrg` directly.

## Two dimensions and irregular graphs

DMRG represents the state along a chain. Interaction range in this ordering and entanglement across its cuts both affect cost. For a two-dimensional or irregular graph, try a bandwidth-reducing ordering:

```python
from qdk_chemistry_alps import interop

order = interop.rcm_ordering(graph)   # bandwidth-reducing permutation
```

The earlier `graph` is already a chain; the complete example below applies this to a shuffled 3×4 lattice, reducing bandwidth from 10 to 4. If reordering an existing calculation, permute every orbital index consistently in the adjacency, integrals and observables. Smaller bandwidth does not guarantee less entanglement or easier convergence; compare solver cost and accuracy.

## Moving models between the two packages

[`examples/interop.py`](../../examples/interop.py) converts a QDK spin operator to ALPS and back, extracts restricted electronic integrals from an ALPS Hubbard model, and demonstrates graph ordering. Run it from the **ALPS repository root**:

```bash
python python/qdk-chemistry-alps/examples/interop.py
```

```text
Heisenberg chain ground state : -9.974309
round trip worst coefficient  : 0.0e+00
extracted integrals           : h1 (4, 4), eri (4, 4, 4, 4)
3x4 lattice bandwidth         : 10 -> 4 after reordering
```

These conversions preserve the supported operator algebra subject to their documented numerical tolerances and coefficient pruning. Integral extraction rebuilds the operator and compares terms; it rejects terms outside the restricted-integral contract, such as pairing or spin-asymmetric interactions. That rejection applies to the integral route, not the broader native Pauli export. Lattice adjacency conversion also does not preserve every geometric attribute. See [Capabilities and limits](../capabilities.md) before assuming a general lossless round trip.

## Next

- [Measure observables](observables.md) on the result
- [Converge a DMRG run](converging-dmrg.md) when there is no exact reference
