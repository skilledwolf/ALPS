# Quantum Wang–Landau SSE

`qwl` samples spin-1/2 isotropic exchange models with the existing Wang–Landau
and multicanonical updates. It uses the shared native Monte Carlo runner,
TOML configuration and ALEA accumulators. On-site terms, anisotropy, sign
problems and self-bonds are rejected because this kernel cannot sample them.
Lattice/model XML libraries remain available through `input.lattice_library`
and `input.model_library`; their parameters belong in `[parameters]`.

```toml
[parameters]
LATTICE = "chain lattice"
MODEL = "spin"
L = 4
J = 1.0
CUTOFF = 24
NUMBER_OF_WANG_LANDAU_STEPS = 8
SWEEPS = 40000
[execution]
seed = 137
chains = 4
[output]
results = "qwl.out.h5"
checkpoint = "qwl.checkpoint.h5"
```

Run `qwl --validate run.toml`, then `qwl run.toml`. `qwl --schema [run.toml]`
prints the application schema, including supplied model/lattice parameters.
The default Zhou–Bhatt refinement and the flatness-based alternative remain
available, as do intermediate coefficient snapshots and optional magnetic
measurements. Omitting `SWEEPS` retains histogram-based production completion.
Nonzero expansion windows start with a valid diagonal operator string.

Independent chains use global chain IDs and seeds in serial or MPI builds.
`execution.rng` selects `mt19937` or `lagged_fibonacci607`. Time limits,
`execution.max_sweeps` and periodic checkpoints use the same runner as
`simplemc`/`spinmc`. Set `input.checkpoint` and distinct output paths to resume;
process counts may change. Increasing `SWEEPS` extends production, including
completed runs, without restarting refinement. Model, window, RNG, chain count
and batch capacity must match. Checkpoints contain the operator string, spins,
refinement/production state, histograms, RNG and exact ALEA accumulator state.

Coefficient/histogram snapshots use mean estimators: one finished run supplies
one estimate, not an independent sample at every visited order. Traversal times
retain native batches and autocorrelation estimates. Each chain's results remain
under `/simulation/realizations/0/clones/<id>/results`; the pooled collection is
under `/simulation/results`. Unvisited orders may have unavailable magnetic or
fraction estimates. The evaluator rejects missing magnetic evidence with nonzero
statistical weight; increase production rather than treating it as zero.

`qwl_evaluate --T_MIN 0.5 --T_MAX 2 --DELTA_T 0.1 qwl.out.h5` generates the
existing XML plot curves. Python's `pyalps.evaluateQWL` reads these plots.
Evaluation uses the final histogram-corrected coefficients for both thermal and
magnetic observables, and evaluates each chain before averaging curves. Turning
off insertion/removal combinatorial factors is accounted for when publishing
coefficients. Positive lower-order windows give conditional energy, heat and
magnetic estimates; they cannot establish absolute free energy or entropy.
All temperature curves retain the finite expansion-window truncation. Check
cutoff convergence when interpreting low-temperature results.

Released QWL final estimates can be migrated offline with
`alps-hdf5-convert old.out.run1.h5 native.h5 --qwl-sites N`, then evaluated with
`qwl_evaluate`. The supplied site count is checked against order-zero
normalization. Pooled logs cannot recover individual chains: use per-run files
when the original count exceeds one. This profile converts statistical results,
not released scheduler checkpoints into resumable native solver state.

The algorithm is described by M. Troyer, S. Wessel and F. Alet,
Phys. Rev. Lett. **90**, 120201 (2003). The original QWL implementation is by
Stefan Wessel; cite it and the ALPS project when publishing results.
