# Native quantum Monte Carlo runs

`worm` and `dirloop_sse` retain their executable names and scientific update
kernels, using the same TOML CLI and native ALEA runner as `qwl` and `spinmc`.
For example:

```toml
[parameters]
LATTICE = "chain lattice"
MODEL = "spin"
L = 4
local_S = 0.5
J = 1.0
T = 1.0
THERMALIZATION = 2000
SWEEPS = 50000
SKIP = 3
[execution]
seed = 137
chains = 3
bins = 64
[output]
results = "run.h5"
checkpoint = "run.checkpoint.h5"
```

Use `worm --validate run.toml`, then `worm run.toml`, or substitute
`dirloop_sse`. `--schema [run.toml]` describes the run fields and supplied model
parameters. Model and lattice XML libraries remain supported through
`input.model_library` and `input.lattice_library`; scheduler XML job files and
scheduler flags are no longer run interfaces. `execution.rng` selects
`mt19937` or `lagged_fibonacci607`.

Set `input.checkpoint` to resume and use distinct output paths. Native checkpoints
retain the physical configuration, RNG, counters, tuning state and unfinished
ALEA batches. Increasing `SWEEPS` extends production, including completed runs.
Other run parameters, chain count and batch capacity must match. Independent chains run in serial or MPI; process counts can change on restart.
The shared runner provides sweep/time limits and periodic checkpoints.

The pooled physical results are under `/simulation/results`; raw per-chain
batches and autocorrelation diagnostics remain under
`/simulation/realizations/0/clones/<id>`. Each signed observable carries an aligned
sign denominator. Nonlinear estimates use weighted joint jackknife analysis;
undefined estimates are omitted with reasons in `/simulation/unavailable`.
Density fluctuations are centered before accumulation to avoid cancellation.
`worm_evaluate run.h5` refreshes derived native results using atomic publication.

`dirloop_sse` retains `minbounce`, `heatbath` and `locopt` directed loops,
local/site compressibility and Green-function measurements, with or without
worm matrix-element weights (`NO_WORMWEIGHT`). The Green function
is the symmetrized equal-time correlation
`(<b†(i)b(j)> + <b(i)b†(j)>)/2` (spin raising/lowering operators for spin models).
Its numerator includes the sign of the open-worm configuration, including
antiferromagnetic correlations in otherwise sign-free systems. Normalization
accounts for `SKIP`. Green measurements combined with nonzero `WORM_ABORT` or
particle/magnetization restrictions reject: those combinations do not have a
correct implemented estimator.

`worm` retains onsite and nonlocal interactions, winding estimators, optional
one-dimensional stiffness and chain densities. `CHAIN_KAPPA` groups sites by
finite transverse coordinates, ordered by coordinate, and normalizes each group
by its own number of sites. Graphs without periodic bond vectors omit winding
and stiffness estimates. Its legacy Green-function implementation was compiled
out; requesting it rejects, as do the unimplemented site-compressibility and
bond-type-stiffness options.

For fixed particle number, `NUMBER_OF_PARTICLES`, positive `CORRECTION` and
`ADJUST` (default `mu`) tune the sampling fugacity. Tuning must finish in the first
quarter of `THERMALIZATION`; increase that budget if necessary. `ADJUST` may only
change onsite energies by a uniform multiple of particle number plus constants,
so relative weights within the selected sector remain unchanged. Energies refer
to the requested Hamiltonian; each chain records its actual
`sampling_parameters`. Production selects the requested sector, including the
vacuum. Conditional number fluctuations, and hence the reported fluctuation
compressibility, vanish in a fixed-number sector.

Released scheduler checkpoints cannot yet be resumed by these executables.
Keep originals and released readers until offline conversion is implemented.
The `loop` migration, released checkpoint converters and tutorial migration remain
part of the ALEA consolidation; the legacy library has not yet been removed.
