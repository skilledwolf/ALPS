# Native classical spin Monte Carlo

`spinmc` reads typed TOML runs, samples with native ALEA batches and writes HDF5
results. It supports Ising, XY, Heisenberg, O(4) and Potts (`q = 3`, `4` or `10`),
with local and legal cluster updates. The private driver shared with `simplemc`
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

Independent `execution.chains` use successive seeds and run serially. The graph
uses one `execution.disorder_seed`, defaulting to the base seed. The RNG is
`mt19937`; old alternate-RNG checkpoints require the released executable.
Multi-process launches and legacy XML run inputs are rejected. Separate TOML
tasks can be distributed by an external job launcher.

Raw batches from all chains are concatenated with their actual sample weights.
Native joint jackknife analysis writes Specific Heat, Connected Susceptibility,
Binder Cumulant U2/U4 and magnetization/Binder slopes directly to the result file;
there is no separate evaluator. Improved cluster estimators have explicit
`Improved ...` names and remain separate from direct moments used in derivatives.
`Connected Susceptibility` uses the radial convention
`beta*N*(<|M|^2> - <|M|>^2)`, where `M` is magnetization per site; it measures
fluctuations of its magnitude. `Susceptibility` is the disconnected second
moment `beta*N*<|M|^2>/d`, with `d` the spin dimension (one for Potts).

Derived results require adequate occupied batches, nonzero denominators and
numerically resolved moment differences. If a difference is within a
conservative floating-point cancellation threshold in the full estimate or any
leave-one-out estimate, the corresponding optional result is omitted. This
prevents very low temperatures from amplifying unresolved cancellation into a
large Specific Heat or slope; raw physical moments and checkpoints remain
available. This screens the final subtraction; it cannot recover precision
already lost inside signed sums. Zero inverse temperature gives exact zero thermal prefactors.
Batch slots are adaptively permuted and result files do not retain chronology.

Python uses `pyalps.run_io.write_run_files` and `execute`, then ordinary
`pyalps.loadMeasurements`; see the [classical tutorials](../../../../tutorials/README.md).
The [format converter](../../../tools/hdf5/README.md) extracts recoverable analysis
evidence from supported released HDF5 profiles. It does not fabricate native
restart state or missing joint covariance.
