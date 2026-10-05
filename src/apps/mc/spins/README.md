# Native classical spin Monte Carlo

`spinmc` reads typed TOML runs, samples with native ALEA batches and writes HDF5
results. It supports Ising, XY, Heisenberg, O(4) and Potts (`q = 3`, `4` or `10`),
with local and legal cluster updates. The SDK runner (`alps::mc`) shared with `simplemc`
handles independent chains, complete checkpoints and pooled analysis.

```toml
[parameters]
MODEL = "Heisenberg"
UPDATE = "cluster"
LATTICE = "chain lattice"
L = 60
T = 1.0
J = [-1.0]
THERMALIZATION = 10000
SWEEPS = 500000

[execution]
seed = 42
chains = 2
bins = 128

[output]
results = "results.h5"
checkpoint = "checkpoint.h5"
```

Run `spinmc --validate run.toml`, then `spinmc run.toml`. Several run filenames
may be passed together; every task and checkpoint is validated before execution.
`spinmc --schema [run.toml]` reports the schema, including graph variables and
numbered couplings found in that run. Paths are relative to the TOML file.
Lattice/graph definitions remain XML scientific resources, supplied through
`input.lattice_library`; run orchestration uses TOML.

Supply exactly one `T > 0` or `beta >= 0` (zero inverse temperature is infinite
temperature). `CONVENTION = "classical"` uses spin length `S` (default one) and
positive `J` favors parallel spins. `CONVENTION = "quantum"` uses
`sqrt(S*(S+1))` (default `S = 0.5`) and reverses the coupling sign; this convention
requires zero field. `S<site type>` overrides site lengths.

`J`, `J<bond type>`, `D` and `D<site type>` are TOML arrays. For an N-component
spin, one value is isotropic, N values specify the diagonal, N*(N+1)/2 values
specify the symmetric lower triangle row by row, and N*N values specify a full
matrix row by row. Onsite `D` must be symmetric. `h` is an array with one value
(along the last axis) or N components; `g` scales the field. Potts requires
isotropic `J` and zero `D`. Its magnetization is color0 occupation, and its field
favors color0: supply `[h]` or `[0, h]`; a nonzero first component is rejected.
Custom scalar lattice bindings retain their types.

`UPDATE = "auto"` chooses cluster updates when the field vanishes, the matrices
are isotropic and the couplings are unfrustrated; otherwise it chooses local
updates. Explicit `cluster` rejects incompatible models instead of sampling
an invalid Hamiltonian. Potts clusters require ferromagnetic couplings.
Thermalization defaults to `SWEEPS/10`; cluster warm-up counts updated sites,
while `SWEEPS` counts subsequent production updates.

`execution.max_sweeps` limits update calls per chain for this invocation,
including warm-up. `execution.time_limit` and `checkpoint_interval` are seconds.
To resume, set `input.checkpoint` to the previous checkpoint and use a distinct
output checkpoint filename. Keep the model, graph, seed, chain count,
batch capacity and thermalization fixed. `SWEEPS` may increase after completion,
but cannot fall below the measurements already taken. Set `THERMALIZATION`
explicitly when changing `SWEEPS` so its default does not change.
Checkpoints retain spins, RNG, progress and unfinished batches; analysis
results cannot replace restart state. `execution.error_variable` and a positive
`error_limit` optionally stop on a scalar native measurement's uncertainty once
at least two batches are occupied. `execution.print_sweeps` sets a diagnostic
interval in updates (zero disables printing). These stopping and printing
options may change when resuming; the current run's options take effect.

Independent `execution.chains` use successive seeds. In an MPI build,
`mpiexec -n 4 spinmc run.toml` distributes their global IDs across ranks; set
`execution.chains` to at least the process count to use every rank. The same
chains have identical streams, raw bins and diagnostics in serial and MPI, and
checkpoints can resume with a different process count. The graph uses one
`execution.disorder_seed`, defaulting to the base seed. RNG choices are
`mt19937` (default) and `lagged_fibonacci607`; both support exact native restart.
Released scheduler seed expansion and buffered checkpoint layouts differ.

The root publishes checkpoints and results. At publication, other ranks transfer
native HDF5 checkpoint bytes through MPI using private temporary spools; this
needs no shared scratch directory and adds no second state schema. Rank failures
propagate at checks every 32 updates. Local budgets/signals are checked each
update. Aggregate checkpoint transport is limited to approximately 2 GiB per
publication by Boost.MPI's count range.

Raw batches from all chains are concatenated with their actual sample weights.
Native joint jackknife analysis writes Specific Heat, Connected Susceptibility,
Binder Cumulant U2/U4 and magnetization/Binder slopes directly to the result file;
there is no separate evaluator. Improved cluster estimators have explicit
`Improved ...` names and remain separate from direct moments used in derivatives.
`Connected Susceptibility` uses the radial convention
`beta*N*(<|M|^2> - <|M|>^2)`, where `M` is magnetization per site; it measures
fluctuations of its magnitude. `Susceptibility` is the disconnected second
moment `beta*N*<|M|^2>/d`, with `d` the spin dimension (one for Potts).

Each batch retains centered physical covariance of `[Energy, |M|, M², M⁴]`.
Specific Heat, Connected Susceptibility and thermal slopes use these moments
with weighted native jackknife, preserving within-bin fluctuations and avoiding
subtraction of large raw energy moments. Constant streams give exact zero
responses, including at very low temperatures. The native covariance accumulators
are saved under each chain's `physical_moments/<bin>` and in its checkpoint.
Undefined denominators, insufficient occupied batches, or unrepresentable
estimates produce reasons under `/simulation/unavailable/<observable>`; raw
evidence remains available. Ordered diagnostics below retain chain chronology.

Python uses `pyalps.run_io.write_run_files` and `execute`, then ordinary
`pyalps.loadMeasurements`; see the [classical tutorials](../../../../tutorials/README.md).
The [format converter](../../../tools/hdf5/README.md) extracts recoverable analysis
evidence from supported released HDF5 profiles. It does not fabricate native
restart state or missing joint covariance.

Each output also retains diagnostics for each independent chain under
`/simulation/realizations/0/clones/<id>`: `autocorrelation/<observable>` holds a
native logarithmic variance hierarchy (kind 4); `series/<observable>` holds
batch sums, weights and chronological offsets in the native accumulator codec.
These measurement states alone are not simulation checkpoints. Sort occupied
bins by offset before plotting; different widths are intentional. Set
`execution.bins` at least as large as `SWEEPS` (and even) to retain every sample
individually, at a proportional memory/disk cost. Histograms of individual
samples require this choice; histograms of coarse bin means describe a different
distribution. Python can use NumPy histogram routines on the retained samples.

`pyalps.loadBinningAnalysis` reads the first chain by default; pass `respath` for
another chain. Its datasets retain `native_result` for level counts, uncertainties
and `tau`. A finite autocorrelation estimate does not prove convergence. Never
join independent chains into one chronological series.
