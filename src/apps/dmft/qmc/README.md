# DMFT run configuration

`dmft --validate run.toml` prepares the numerical input and checks the selected
solver's settings without launching it or writing results. `dmft run.toml` runs
the self-consistency loop. Run files contain values; the application supplies
its schema, identity and version. Legacy parameter files and parameter archives
are no longer accepted as run configurations.

```toml
[parameters]
BETA = 2.0
U = 0.0
MU = 0.0
t = 1.0
N = 16
NMATSUBARA = 8
SWEEPS = 3000
THERMALIZATION = 100
N_MEAS = 4

[output]
results = "results.h5"

[execution]
solver = "hybridization"
loop = "omega"
max_iterations = 2
seed = 42
```

`execution.solver` selects the `hybridization` (CT-HYB), `interaction` (CT-INT)
or `hirschfye` executable in the ALPS `bin` directory (`ALPS_BIN_PATH`). Any other
value names a custom executable, resolved against that directory unless it is
absolute. The driver
maps its time/frequency grids and flavors to the child application's parameters.
It passes one TOML file and a separate temporary HDF5 archive holding one vector
per flavor, `/Delta_<f>` for hybridization solvers and `/G0_<f>` otherwise;
temporary files are removed after success or failure. Each child process uses
the requested seed. `execution.time_limit = 0` disables the child time limit;
finite `SWEEPS` still terminate sampling.
DMFT launches these solvers serially. Direct standalone MPI runs use independent
chains, distinct seeds and aggregate `SWEEPS`.

`dmft --schema run.toml` prints the composed application schema, including the
selected solver's scientific settings. The embedded DMFT base schema is extended
by the CT-HYB, CT-INT or Hirsch-Fye schema and the selected flavor count. A
custom solver requires `input.solver_schema` and `execution.solver_input`, either
`"delta"` for a hybridization solver or `"g0"`, and the same one-file protocol.
The tau loop requires a hybridization solver; the others use the omega loop.

Paths resolve relative to the run file. Optional `input.initial_omega` and
`input.initial_tau` supply the existing numerical Green-function text layouts:
a coordinate followed by each flavor's value, with complex values written as
`(real,imag)`. Dimensions and finite values are checked. `input.dos` supplies
energy/density pairs for each band on an increasing, uniform, odd grid of at
least three points. It selects the general omega Hilbert transform; the previous
`SEMICIRCLE_HILBERT = false` default is preserved. `TWODBS`, `tprime` and `L`
retain the existing two-dimensional bandstructure settings.

`input.interaction_matrix` retains the sparse DMFT text convention `i j U_ij`.
CT-INT accepts a finite, symmetric matrix with zero diagonal, including isolated
flavors or an entirely zero matrix. Without this file, `U`, `J` and optional
`"U'"` (default `U - 2*J`) assemble the interactions for paired flavors.
For segment CT-HYB, optional `input.retarded_interaction` and its format/coordinate
settings follow the [CT-HYB scientific-input contract](hybridization/README.md).
These files contain scientific data, not run settings.

Results retain `/simulation/iteration/<n>/results` and `/simulation/results`.
Typed scientific parameters are stored at `/parameters`; `/run_config` records
resolved input/output/execution settings and input/default/derived provenance.
`pyalps.loadDMFTIterations` can analyse these results. A saved configuration is
not a complete solver restart checkpoint.

Standalone `interaction` uses modern ALEA batches. Physical measurements retain
their signed numerator and sign together, so W and density estimates are
normalized by average sign with their covariance retained. MPI collection pools
the joint batches before evaluating ratios, including partial bins and empty
ranks. `execution.bins` is the even number of local batch slots, at least two;
it does not discard measurements or repartition collected replicas. A single
occupied bin has an unavailable error. A zero sign denominator, including a
singular jackknife leave-out estimate, rejects the analysis.

Its `/simulation/results/<name>` groups use the versioned native ALEA batch
codec (`@version=1`, `@kind=5`). `pyalps.loadMeasurements` retains the reported
means and errors; `pyalps.alea.BatchResult.read` also exposes per-bin sums and
counts. Green-function paths are unchanged. The complete HDF5 output replaces
the destination only after successful evaluation and serialization. These
analysis files do not resume the solver's Markov chain.

CT-INT supports between 2 and 128 flavors with one density-density kernel. It proposes
uniformly among unordered interacting flavor pairs, preserving detailed balance
for sparse matrices without retry loops on zero rows. Its joint signed batches
retain all `FLAVORS * FLAVORS` density products, including `n_f² = n_f`.
General interaction-dependent Fourier tails use those density moments and each
flavor's full interaction row. `PertOrder` supplies order diagnostics in HDF5.
`interaction --schema run.toml` and Python's `ctint.schema(parameters={...})`
include the selected flavor count's `EPS_<f>` and `EPSSQ_<f>` settings.

Standalone CT-INT's `input.atomic = true` uses a zero-energy bare level. It
derives zero `EPS_<f>` and `EPSSQ_<f>` Fourier moments; explicit nonzero moments
are rejected. This keeps both measurement paths consistent with the atomic input.

Text sidecars require `output.text = true` and an existing
`output.text_directory`; final numerical Green files require explicit
`output.final_tau` or `output.final_omega` paths. Outputs cannot replace run or
scientific input files. The driver neither changes the working directory nor
rewrites the user's run file.

Python callers write run files with `pyalps.run_io.write_run_file` and run them
with `pyalps.run_io.execute("dmft", run_files)`, which also accepts job manifests.
All runs are validated before execution. A successful call returns the absolute
`output.results` path of each run; process failures raise an exception.

Standalone `hirschfye` has its own schema (`hirschfye --schema`, installed under
`share/alps/schemas`) and reads `/G0_<f>` vectors from the HDF5 file `input.g0`.
DMFT selects this executable with `execution.solver = "hirschfye"`.

Hirsch-Fye uses native ALEA batches for Sign and the two joint Green/sign
streams. Every physical time sample includes the beta endpoint before
accumulation; ratios, means and errors therefore follow one convention.
Replicas pool raw measurements before weighted jackknife analysis, preserving
partial batches. `execution.bins` must be even and at least two (default 128).
Warm-up sweeps are excluded, including the transition sweep. Native kind-5
results, `/G_tau`, `/G_omega`, typed parameters and provenance are published
together through checked atomic HDF5 replacement. Python's ordinary
`loadMeasurements` reads the means and errors. These files do not restart the
auxiliary-spin chain. The disabled four-point implementation and its
`MEASURE_FOURPOINT_FUNCTION` / `FOURPOINT_INTERVAL` settings are removed.
