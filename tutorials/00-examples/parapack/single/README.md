# Native heat-bath Ising and temperature scans

`ising_single` uses the shared native MC runner, consolidated parameters and
TOML. Build it from the parent examples project against an installed ALPS SDK:

```sh
build/parapack/single/ising_single --schema
build/parapack/single/ising_single --validate parapack/single/run.toml
build/parapack/single/ising_single parapack/single/run.toml
build/parapack/single/ising_single parapack/single/scan.toml
```

Specify exactly one `LATTICE` or `GRAPH`; an optional `input.lattice_library`
selects another XML library. The Hamiltonian is `-J*sum_b s_i*s_j`, with spins
±1 and uniform `J` (either sign). Each sweep updates all graph-color classes
with the heat-bath flip probability. Neighboring distinct sites have different
colors. Self-loop bonds contribute a constant energy and do not enter a flip's
energy change. The logistic probability is evaluated without exponential
overflow and remains valid for uncoupled spins at very large inverse temperature.

With `ALGORITHM="ising; temperature scan"`, the temperatures are
`INITIAL_TEMPERATURE + i*DIFF_TEMPERATURE`, for `i=0,...,NUM_TEMPERATURES-1`.
The initial stage uses `INITIAL_THERMALIZATION` when supplied; later stages use
`THERMALIZATION`. Production `SWEEPS` are measured separately at each temperature.
The **same physical configuration and RNG stream continue between stages**.
Fixed-temperature runs use `ALGORITHM="ising"` (the default) and require `T`.
All temperatures and inverse temperatures must be finite and positive.

`execution.chains` selects independent configurations or scans; an MPI-enabled
SDK distributes these through the shared runner. `execution.seed` supplies the
base seed, and `execution.rng` selects MT19937 or lagged Fibonacci.
`execution.max_sweeps` budgets sweeps per invocation and chain, including all
warmups. Native checkpoints preserve spins, ordered topology, RNG state, the
completed-sweep count, and every stage's native accumulators and unfinished bins.
Continue with `input.checkpoint` and new output paths. A fixed-temperature run
can extend `SWEEPS` with explicit unchanged warmup. A scan's stage lengths must
match, since changing them would move the temperature-transition boundaries.

Fixed-temperature results use `/simulation/results`. Scans reuse the native
ensemble layout `/simulation/replicas/<stage>/results`, with each stage's `T`
in its own parameters. `pyalps.loadMeasurements` and `loadBinningAnalysis` load
these datasets directly. Native joint samples retain the original extensive
energy and magnetization moments and number of sites; derived output includes
heat capacity per site and the original Binder ratio `<M^2>^2/<M^4>`.
Heat capacity uses centered physical moments within the same bins to avoid
subtracting large nearly equal raw energies. All independent chains and partial
bins participate in weighted jackknife analysis. Each chain also retains a
vector autocorrelation hierarchy, ordered as site count, energy, squared energy,
magnetization, squared magnetization and fourth magnetization moment.
Insufficient statistics or undefined estimates are recorded explicitly.

Optional OpenMP updates within each color class remain supported. Enable
`ALPS_ENABLE_OPENMP` when configuring the examples (using an SDK's compatible
OpenMP runtime), and set `OMP_NUM_THREADS`. Random proposals are generated in a
fixed serial order; only independent spin updates run in parallel. Results and
checkpoints are therefore identical across thread counts. Small color classes
execute serially to avoid thread-launch overhead.

`kernel.hpp` also supplies the physical update to the temporary legacy adapter
in `ising.h`, still used by the unported exchange example. The native spatial
command reuses this model's statistics and scan handling with a distributed kernel.
The native command does not include that adapter. Finish runs from released Parapack checkpoints with ALPS 3.0. The offline
converter migrates results for analysis, not physical restart state.
