# Native simplemc

`simplemc` runs the Ising, XY and Heisenberg local-update algorithms on ALPS
lattices or explicit graphs. It reads typed TOML runs and uses native ALEA
batches for statistics and HDF5 checkpoints. A single implementation shares the
graph, coupling, measurement, checkpoint and snapshot code across the models.

```toml
[parameters]
ALGORITHM = "ising"
LATTICE = "square lattice"
L = 16
T = 2.3
J = 1.0
H = 0.0
SWEEPS = 65536
THERMALIZATION = 8192

[output]
results = "ising.h5"
checkpoint = "ising-checkpoint.h5"
snapshot_prefix = "spins"

[execution]
seed = 42
chains = 2
bins = 128
snapshot_interval = 16384
```

Run `simplemc --validate run.toml` to check a run without writing files, then
`simplemc run.toml` to execute it. `simplemc first.toml second.toml` validates
every task and its output paths before starting any task. `simplemc --schema
run.toml` prints that run's schema, including its scalar graph bindings.

Specify exactly one `parameters.LATTICE` or `parameters.GRAPH`. Built-in lattice
XML is found through the installed resources or `ALPS_XML_PATH`; custom XML goes
in `input.lattice_library`. Scalar lattice bindings such as `L`, `W`, `a` and
custom extent/coordinate variables retain their native TOML types. `J0`, `J1`,
and other `J<bond type>` keys override `J` for those bond types. The XY field is
along X; the Heisenberg field is along Z. Omit `T` for infinite temperature.
Omitted thermalization defaults to `SWEEPS / 8`.

Each of `execution.chains` independently completes `THERMALIZATION + SWEEPS`
sweeps, using seed `execution.seed + chain index`. These are chains of the same
Hamiltonian and graph. `execution.disorder_seed` sets the graph's quenched
disorder/depletion seed and otherwise defaults to `execution.seed`; all chains
share it. Use separate explicit tasks for different disorder realizations. The
driver distributes global chain IDs across ranks in MPI builds; for example,
`mpiexec -n 4 simplemc run.toml`. Set `execution.chains` to at least the process
count to use every rank. Streams and raw evidence are identical in serial/MPI,
and checkpoints can resume with a different process count. The supported
RNG choices are `mt19937` (default) and `lagged_fibonacci607`. Both retain
exact native checkpoint continuation. Seeding uses Boost’s integer seed constructor;
released scheduler seed expansion and buffered checkpoint layouts differ.

Raw native batches are concatenated across chains with their sample weights.
Specific Heat uses centered energy moments retained within every batch, avoiding
subtraction of large raw moments and preserving physical fluctuations within
bins. Native weighted jackknife supplies its uncertainty. Constant energy gives
exactly zero heat capacity, including at very low temperatures. Binder ratios
use aligned magnetization batches: `<m²>² / <m⁴>`, including X/Z projections.
Undefined denominators, insufficient occupied batches, or unrepresentable
estimates produce explicit reasons under `/simulation/unavailable/<observable>`.
Raw evidence remains available. Each chain's `physical_moments/<bin>` contains
the native covariance accumulator needed to reproduce thermodynamic analysis;
these moments also resume exactly from the application checkpoint.

`execution.time_limit` stops after that many seconds, and
`execution.max_sweeps` supplies a per-chain sweep budget for this invocation;
zero means run to completion. Both include thermalization. Signals also stop
the run. `output.checkpoint`, when supplied, is written atomically at stopping
and periodically according to `execution.checkpoint_interval` (seconds; zero
disables periodic saves). Results from the completed measurements are published
even when interrupted during warm-up.

Resume by setting `input.checkpoint` to the saved file, keeping the model,
graph, seed, chain count, bin count and thermalization unchanged. `SWEEPS` may
increase, including after completion, but cannot fall below the measurements
already taken. Set `THERMALIZATION` explicitly when changing `SWEEPS` so its
default does not change. Stopping budgets and output paths may change.
Use a distinct output checkpoint filename. Checkpoints
contain spins, total sweep count, RNG and live native accumulators, including
partially filled batches. All chains are loaded and checked before any output
is written.

Snapshots are direct ASCII VTK files named
`<snapshot_prefix>.clone<chain>.<total sweep>.vtk`, with chain indices starting
at one. Prefix namespaces cannot overlap task inputs or outputs. Existing
snapshots are preserved; use a fresh prefix when a restart would repeat them.
The [snapshot tutorial](../../../../tutorials/03-mc/09-snapshot/) generates files
compatible with its ParaView state files. `snap2vtk` remains an offline converter
for snapshots written in the released XDR format.

The old XML front end and Parapack worker/checkpoint formats are removed from
this executable. Old `NUM_CLONES`, `SEED`, `RNG`, `SNAPSHOT_INTERVAL` and
`LATTICE_LIBRARY` settings move to the typed fields above. Evaluate legacy
parameter expressions before writing TOML; arrays represent values, and task
sweeps are explicit run files. Old statistical HDF5 bins can be converted with
the [offline converter](../../../tools/hdf5/README.md); that analysis conversion
does not reconstruct a resumable simulation checkpoint.

Python orchestration uses the same executable:

```python
from pyalps.run_io import execute, write_run_files

job = write_run_files("temperatures", [
    {"parameters": {"ALGORITHM": "ising", "LATTICE": "chain lattice",
                    "L": 16, "T": temperature, "SWEEPS": 4096},
     "output": {"results": f"temperature{i}.h5"}}
    for i, temperature in enumerate((1.5, 2.0, 2.5))
], baseseed=42)
results = execute("simplemc", job)
```

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

MPI publishes on the root after transporting native HDF5 checkpoint bytes using
private temporary spools; ranks need no shared scratch directory. Failures reach
every rank at checks every 32 updates, with local stopping checked every update.
Boost.MPI limits aggregate checkpoint transport to approximately 2 GiB per
publication. The same scientific checkpoint codec serves serial and MPI runs.
