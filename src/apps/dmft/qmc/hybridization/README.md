# Segment CT-HYB runs

The executable reads a TOML run file. `hybridization --validate run.toml`
checks the configuration and numerical input without sampling or writing results;
`hybridization run.toml` runs it. The application supplies its identity and schema
version. The installed `share/alps/schemas/cthyb.toml` describes supported keys,
types, defaults, and bounds.

```toml
[parameters]
BETA = 10.0
N_ORBITALS = 2
N_TAU = 100
N_MEAS = 10
THERMALIZATION = 1000
SWEEPS = 10000
U = 2.0
MU = 1.0

[input]
delta = "delta.dat"

[output]
results = "result.h5"

[execution]
seed = 42
```

Paths in a file resolve relative to that file, independently of the working
directory. `execution.time_limit = 0`, the default, disables the time limit;
finite `SWEEPS` still terminate the run. Positive limits are in seconds.
Text results require `output.text = true`; `output.text_directory` selects an
existing destination directory. `output.base_path` selects an HDF5 group for
scientific results. The archive also saves typed scientific parameters and the
resolved four-section configuration, including provenance. HDF5 publication
is atomic: a failed save preserves the previous archive. Explicit text files
are separate outputs and do not participate in that transaction.

`execution.bins` selects an even number of native ALEA batch slots, at least
two (default 128). It replaces the scientific `NUM_BINS` parameter.
Sign-weighted time measurements use the sign averaged over `N_MEAS` updates;
frequency, Legendre and other final-configuration measurements use that
configuration's sign. Independent MPI chains pool raw measurements before
normalization, preserving partial batches. `SWEEPS` remains aggregate work
across MPI ranks; each rank uses its own seeded stream.

Sign-weighted bounded measurements retain joint numerator/sign batches and
weighted jackknife errors. Sign and order diagnostics use ordinary batches.
Large G2/H2 tensors retain componentwise numerator/sign
covariance in linear memory; their errors assume independent samples and do
not include autocorrelation. Results under `simulation/results` use the
canonical ALEA batch (kind 5) or variance (kind 2) codec, readable by
`pyalps.loadMeasurements`. G/F time endpoints, errors and covariance follow
one convention; covariance is always that of the mean and the redundant
`ACCURATE_COVARIANCE` selector is removed.

Text Delta input has `N_TAU + 1` rows, with a coordinate and one nonpositive value
per orbital. Coordinates are consecutive indices starting at zero by default.
Set `input.delta_coordinate = "tau"` for a uniform grid from zero to `BETA`.
For HDF5, set `input.delta_format = "hdf5"` and supply dimension-checked vectors
`/Delta_0`, `/Delta_1`, etc.

Optional `input.interaction_matrix`, `input.chemical_potential`, and
`input.retarded_interaction` replace the corresponding scalar model input or add
retardation. Their format keys are `interaction_format`,
`chemical_potential_format`, and `retarded_interaction_format`, respectively.
Text matrix/vector input contains exactly the required number of values. HDF5
datasets remain `/Umatrix`, `/MUvector`, and `/Ret_int_K` plus `/Ret_int_Kp`.
Text retarded input has a coordinate, K, and K-prime per row, with an optional
`retarded_interaction_coordinate = "tau"`. Numerical files remain supported;
HDF5 archives of legacy run parameters are no longer accepted as run input.

Native callers use `alps::cthyb::prepare_run(run)` and
`alps::solvers::cthyb(run)`. Python callers use `pyalps.cthyb.prepare` with
separate `parameters`, `input`, `output`, and `execution` dictionaries, then
`pyalps.cthyb.solve(run)`; `pyalps.run_config.load(file, pyalps.cthyb.schema())`
loads a TOML run file for `solve`.
The current solver does not serialize its segment configuration for restart;
result files do not provide a restart checkpoint.
