# C++ Ising simulation

Implement `step()` and `observables()` in `ising-skeleton.cpp`: random-site
Metropolis updates on a periodic square lattice, then the per-site energy and
magnetization moments. The completed physics is in `solution/ising.cpp`.
Both use `simulation.hpp` for typed parameters, independent RNG state, native
ALEA statistics and atomic HDF5 output.

With an installed ALPS SDK, configure and build the solution from this directory:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/path/to/alps/install
cmake --build build --config Release
```

After completing the exercise, build it with
`cmake --build build --config Release --target ising-exercise`.
The executable names remain `ising` and `ising-exercise`. With no arguments,
either runs the original L=16, beta=0 through 1 scan (2500 warmup and 5000
production sweeps per point) in the current working directory.

For a specific run, save this as `run.toml`:

```toml
[parameters]
L = 16
BETA = 0.4
THERMALIZATION = 2500
SWEEPS = 5000
[execution]
seed = 42
bins = 64
[output]
results = "ising.h5"
```

Run `build/ising --validate run.toml` to validate without sampling or writing,
then `build/ising run.toml`. `--schema` prints the parameter schema. Multiple
run files are supported; their output paths must be distinct. Relative paths
are resolved from each TOML file. `bins` must be even and at least two.

The five original observables (`E`, `m`, `|m|`, `m^2`, `m^4`) are native batch
results under `/simulation/results`. Their aligned samples are also retained
in `/simulation/joint`, so subsequent nonlinear analysis can use their mutual
correlations. Native per-observable autocorrelation levels are stored under
`/simulation/realizations/0/clones/0/autocorrelation`. Use
`pyalps.loadMeasurements` and `pyalps.loadBinningAnalysis` to inspect them.
Printed errors come from batches; the convergence heuristic and available tau
come from the autocorrelation hierarchy. Short runs need not establish an
error plateau. These files contain analysis results, not resumable spin states.

Released tutorial files can be migrated offline with `alps-hdf5-convert`
(`--alea-batches`, `--alea-summary`, or `--alea-autocorr` per observable,
depending on the retained evidence); see the converter documentation.
