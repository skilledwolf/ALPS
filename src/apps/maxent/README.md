# MaxEnt run configuration

Run `maxent --validate run.toml` to load data, apply defaults and check input without
creating output. Run `maxent run.toml` to calculate the spectrum. Legacy parameter
files and HDF5 archives are no longer accepted as run configurations.

```toml
format_version = 1
application = "maxent"
schema_version = 1

[parameters]
BETA = 2.0
NFREQ = 100
OMEGA_MAX = 4.0

[input]
values = [-0.5, -0.3, -0.3, -0.5]
errors = [0.01, 0.01, 0.01, 0.01]

[output]
results = "spectrum.h5"
text = false

[execution]
time_limit = 60
```

The application owns [schema/maxent.toml](schema/maxent.toml), which is embedded
at build time and installed under `share/alps/schemas`. The run file contains
values, not validation rules. Unknown keys and unsafe conversions are errors.
Relative paths are resolved against the run file's directory. Scientific
cross-field checks (such as consistent `T`/`BETA` or covariance shape) remain
native application code; they apply equally to C++ and Python callers.

Instead of inline arrays, use `[input] data="measurements.h5"` with
`format="hdf5"`. Dataset defaults are `/Data` and `/Error`; override them with
`values_dataset` and `errors_dataset`. Optional `covariance_dataset` and
`tau_dataset` select data in the same file. Covariance is a flat, row-major square
array matching the measurement count. For `format="text"`, each data record is
`index value error`, with indices starting at zero. These files contain scientific
data only. A separate `covariance_file` contains `row column value` records.

`DEFAULT_MODEL="tabulated"` uses `prior_omega` and `prior_density` arrays or a
`prior` file containing `omega density` records. The coordinates must increase
and cover the frequency range. Analytic priors retain their original names and
parameters, listed in the schema.

C++ callers include `<alps/maxent.hpp>` and call
`alps::solvers::maxent(parameters, data, output_file, time_limit, text_output)`;
the return value indicates whether a result was produced. Python callers use
`pyalps.cxx.maxent_c.AnalyticContinuation(parameters, input, output_file,
time_limit=60, text_output=False)`. Both paths copy and validate input before
constructing the numerical solver. Numerical algorithms are unchanged.

Results retain their spectrum datasets and add `/spectrum/widths`. Resolved
parameters use the `alps.params.v1` checkpoint format at `/parameters`; the CLI
also saves `/run_config` with input/output/execution settings and each value's
origin (`input`, `default` or `derived`). Cancellation before calculation leaves
existing results untouched. Text sidecars are opt-in and use the output file's
stem as their prefix.
