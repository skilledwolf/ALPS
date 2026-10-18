<!--
Copyright (c) 2026 Tobias Wolf
SPDX-License-Identifier: MIT

Adapted from skilledwolf/alpsx (private repository), branch
agent/v3-execution-and-api-audit at d37e51227,
python/qdk-chemistry-alps/docs/reference/settings.md. Changes: the settings of
alps_classic_dmrg (kept states, warm-up, fixed-state sweeps) replace those of
alps_dmrg (seed, MPO tolerance, mps container); alps_ed has no engine selector.
-->

# Settings

Configure a calculator before calling `run()`. QDK locks its settings when the run starts; create a new calculator to use different settings in a later calculation.

```python
import qdk_chemistry_alps  # registers the calculators
from qdk_chemistry.algorithms import create

calc = create("multi_configuration_calculator", "alps_classic_dmrg")
calc.settings().set("max_bond_dim", 256)
```

## Shared

| Setting | Type | Default | Effect |
| --- | --- | --- | --- |
| `integral_cutoff` | float | `1e-12` | Integrals below this magnitude are dropped when building the operator (before the symmetry check, in `alps_classic_dmrg`). Raising it shortens the term list at the cost of accuracy. |
| `ci_coefficient_cutoff` | float | `0.0` | Determinants with coefficient magnitude below this threshold are omitted from the returned vector, which is renormalized. `0.0` keeps all of them. Affects only the returned wavefunction, never the energy. |

A positive CI cutoff changes the returned state: its observables need not reproduce the solver's energy or measurements, and `alps_classic_dmrg` then attaches no RDMs to the container. Keep `0.0` for comparisons that need all amplitudes.

`alps_ed` has no other settings: it diagonalizes in the $(N, S_z)$ sector of the electron counts with `pyalps.ed`.

## `alps_classic_dmrg`

| Setting | Type | Default | Effect |
| --- | --- | --- | --- |
| `max_bond_dim` | int | `200` | Kept states in every half sweep, and the cap under error-based truncation. See [Converge a DMRG run](../how-to/converging-dmrg.md). |
| `warmup_states` | int | `20` | Kept states while the blocks grow. |
| `sweeps` | int | `20` | Sweeps before the final measurement pass. Raise for strongly correlated or two-dimensional systems. |
| `convergence_tol` | float | `1e-10` | Lanczos tolerance of each step; not an early-stopping threshold for the run. |
| `truncation_tol` | float | `0.0` | Zero keeps the fixed-state schedule; a positive value targets that discarded weight, capped at `max_bond_dim`. |
| `wavefunction_container` | str | `"auto"` | `"auto"`: a `StateVectorContainer` while the CI vector fits `ci_memory_limit_gb`, otherwise an `RdmContainer`. `"statevector"`: always the CI vector (`MemoryError` beyond the limit). `"rdm"`: an `RdmContainer` without a CI vector. See [Work past the CI limit](../how-to/large-active-spaces.md#choosing-explicitly). |
| `ci_memory_limit_gb` | float | `4.0` | Budget for the CI vector, estimated at 1 kB per determinant of the sector, so about four million determinants. The choice is made before the solve. |

The calculator takes restricted real integrals; it rejects integrals that break the eightfold symmetry.

## Reading settings and results back

```python
value = calc.settings().get("max_bond_dim")
```

After a run, `calc.last_result` holds the `pyalps.classic_dmrg` result: `energy` (with the core energy), `iterations` (direction, energy, truncation error and entropy of each step), `final_truncation_error`, the RDMs (`one_rdm` and `two_rdm` spin-traced, `one_rdm_alpha`, `one_rdm_beta`, `two_rdm_aa`, `two_rdm_ab` and `two_rdm_bb` by spin block) and `observables`. `calc.observables`, set before `run()`, names further operators to measure in the final sweep; see [Measure observables](../how-to/observables.md#what-you-can-measure).

Reaching the configured `max_bond_dim` is a reason to test a larger value, not by itself a measure of the error.

An `RdmContainer` returned by `alps_classic_dmrg` records the run in its metadata, which is serialized with the container:

```python
metadata = wavefunction.get_container().get_metadata()
print(metadata["energy"], metadata["max_bond_dim"], metadata["sweeps"], metadata["final_truncation_error"])
```
