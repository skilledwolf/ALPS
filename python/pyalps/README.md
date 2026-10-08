# pyalps

Python applications and libraries for the Algorithms and Libraries for Physics Simulations (ALPS) project. Binary wheels are available from PyPI:

```sh
python -m pip install pyalps
```

Matplotlib plotting helpers are included with `pyalps`. Install `pyalps[mpi]` for the mpi4py-backed `pyalps.mpi` compatibility layer.

Wheels bundle simulation applications for `pyalps.run_io.execute`, but do not add them to your shell's `PATH` or include auxiliary tools such as `printgraph` and `alps-xml`. Install the [C++ SDK and tools](../../CONTRIBUTING.md#build) for those command-line workflows. Tutorial files are available in the [source collection](../../tutorials/README.md), with an optional SDK installation component.

The bindings are built as a standalone `scikit-build-core` project using nanobind. A source build requires GIL-enabled CPython 3.11 or newer, CMake 3.27 or newer, Ninja, a C++17 compiler, BLAS/LAPACK, HDF5, and an installed shared ALPS C++ SDK. The SDK's numeric version must match `ALPS_VERSION.txt`; CMake rejects a mismatch before compiling. Point `ALPS_DIR` at the SDK's `share/alps` package directory.

Before building the SDK below, follow the [CMake and Ninja setup](../../CONTRIBUTING.md#install-cmake-and-ninja) if your system CMake is older than 3.27.

The `distribution` CMake preset builds the SDK exactly as the wheel CI does. From the repository root:

```sh
cmake --preset distribution
cmake --build --preset distribution

export ALPS_DIR="$PWD/_build/distribution/install/share/alps"
python -m pip install build
python -m build --wheel python/pyalps
```

The wheel is written to `python/pyalps/dist` and can be installed with `python -m pip install`. With ccache installed, configure with `cmake --preset distribution -DCMAKE_CXX_COMPILER_LAUNCHER=ccache` and set `CMAKE_ARGS="-DCMAKE_CXX_COMPILER_LAUNCHER=ccache"` for the wheel build to speed up rebuilds.

`PYALPS_BUILD_SOLVERS=ON` is the default and builds the MaxEnt, CT-HYB, and CT-INT extension modules. Set it to `OFF` through CMake configuration for a smaller core-only developer build. These modules link `ALPS::maxent`, `ALPS::cthyb`, and `ALPS::ctint` from an SDK built with `ALPS_BUILD_APPLICATIONS=ON`; they do not compile solver implementations.

All binding sources and build helpers live under this package directory. The source distribution includes shared release metadata from `cmake/ALPS_VERSION.txt` and `LICENSE.txt`, stored at its root as `ALPS_VERSION.txt` and `LICENSE.txt`. It can be built outside the checkout against an installed SDK, without the C++ source tree.

The [Ising extension example](examples/ising/README.md) demonstrates exporting a simulation through nanobind against the installed SDK and pyalps wheel.

`PYALPS_BUNDLE_APPLICATIONS=ON` independently controls inclusion of the ALPS command-line programs (`spinmc`, `dmrg`, `sparsediag`, `loop`, `qwl`, ...) in `pyalps/bin`. Required shared libraries are included even when this option is `OFF`. `pyalps.run_io.execute` first honors programs already on `PATH`, then search `ALPS_BIN_PATH` and the bundled directory. A bindings-only installation can therefore use separately installed programs through either environment variable. XML resources come from the package unless `ALPS_XML_PATH` is set; installations contain no fallback paths to the machine that built the wheel.

## Editable development

Build and install the C++ SDK using [CMake](../../CONTRIBUTING.md#build), then use the Python environment of your choice for an editable installation. Keep the SDK and bindings on the same compiler, architecture and native dependency stack.

After installing the matching SDK, install build dependencies in your active Python environment and use pip's editable mode:

```sh
export ALPS_DIR="$PWD/_build/distribution/install/share/alps"
python -m pip install "scikit-build-core>=1.0" "nanobind==2.15.0" "cmake>=3.27" ninja \
  "patchelf>=0.14; sys_platform == 'linux'"
python -m pip install --no-build-isolation -e python/pyalps \
  --config-setting "build-dir=$PWD/_build/manual-python"
```

The example uses the SDK installed by the `distribution` preset above. For another SDK, change `ALPS_DIR` to its installed `share/alps` directory. The Linux `patchelf` dependency is required when bundling applications; `--no-build-isolation` makes installing build dependencies your responsibility. Pip installs NumPy, h5py, SciPy, and Matplotlib as runtime dependencies.

For a smaller SDK and core-only editable Python installation on Linux/macOS, use the `sdk` preset and disable both solver bindings and bundled programs. After installing the build dependencies above:

```sh
cmake --preset sdk
cmake --build --preset sdk --parallel 2
cmake --install _build/sdk
export ALPS_DIR="$PWD/_build/sdk/install/share/alps"
python -m pip install --no-build-isolation -e python/pyalps \
  --config-setting "build-dir=$PWD/_build/manual-python-core" \
  --config-setting cmake.define.PYALPS_BUILD_SOLVERS=OFF \
  --config-setting cmake.define.PYALPS_BUNDLE_APPLICATIONS=OFF
```

Python edits take effect in a new interpreter without reinstalling. Rebuild and reinstall the SDK after SDK C++ changes, then rerun the editable install after changing binding C++ sources, build configuration, or packaged runtime resources. If you switch SDKs, compilers, or dependency providers while reusing a manual binding build directory, add `--config-setting cmake.args=--fresh` once to clear cached discovery. Scikit-build-core owns the Python source mapping; CMake installs the native modules, libraries, XML, downstream headers, and CMake package. Resource lookup uses the installed runtime manifest, so it also works when these files live separately from the Python sources.

## Native runtime layout

The extensions share one nanobind library, named `pyalps_nanobind` to avoid filename collisions with other packages, and one packaged copy of each ALPS runtime component (`alps`, `alps_params`, `alps_run_config`, `alps_hdf5`, `alps_utilities`, `alps_osiris`, `alps_xml` and `alps_cli`). The SDK lists these targets in `ALPS_RUNTIME_TARGETS`. On Unix, CMake installs relative runtime paths with the package's `lib` directory first and derives any additional dependency directories from resolved link targets. On macOS, manifest generation also redirects dependencies between the packaged ALPS libraries, so they do not load a second SDK copy. Linux and macOS wheels intended for redistribution must then be repaired with auditwheel or delocate, respectively, to bundle external dependencies and replace build-machine paths; the wheel CI performs this step. See [downstream native extensions](#downstream-native-extensions) for manifest finalization after repair.

## Python compatibility

pyalps supports GIL-enabled CPython 3.11 and newer. Wheels use the ABI of the interpreter that built them; they do not use the limited API or `abi3`. Free-threaded Python is unsupported.

Build dependencies pin nanobind 2.15.0. The installed `runtime.json` records both its version and internals ABI. The downstream CMake package checks that ABI before creating `pyalps::runtime`, so an incompatible nanobind installation produces a configuration error rather than an interpreter abort. Rebuild native consumers against the same SDK and nanobind release as the installed pyalps package.

## Versioning

pyalps does not carry a version of its own. The numeric version is read from `cmake/ALPS_VERSION.txt` — the same file `cmake/ALPSVersion.cmake` reads for `ALPS_VERSION_CORE` — so a release bump is one edit rather than two that can drift. `tests/pyalps/test_wheel_payload.py` fails if the installed version and that file disagree.

A prerelease label cannot live in that file: `project(VERSION ...)` rejects a non-numeric version, and neither `find_package()` matching nor the library SOVERSION has a notion of prerelease ordering. CMake takes it from the `ALPS_VERSION_PRERELEASE` cache variable; the Python build takes it from the environment variable of the same name, using the same vocabulary:

| `ALPS_VERSION_PRERELEASE` | version with `ALPS_VERSION.txt` = 2.3.4 |
|---|---|
| unset | `2.3.4` |
| `beta.1` | `2.3.4b1` |
| `alpha.2` | `2.3.4a2` |
| `rc.1` | `2.3.4rc1` |
| `dev.3` | `2.3.4.dev3` |

In GitHub release builds, the provider takes the prerelease label from `GITHUB_REF` (for example, `refs/tags/v3.0.0-beta.1`). It rejects a tag whose numeric version differs from `ALPS_VERSION.txt`, or whose label conflicts with an explicit `ALPS_VERSION_PRERELEASE`. Wheels and source distributions use the same provider. `python python/pyalps/_build_support/alps_version.py` prints the version a build would produce, from any working directory. An sdist preserves its recorded version when rebuilt without the original build environment.

Note the consequence: because the number is inherited, a Python-only API change cannot be signalled in the pyalps version alone — it takes a bump of `cmake/ALPS_VERSION.txt`, which moves the whole project.
## Downstream native extensions

The C++ SDK supplies `ALPS::alps` for standalone programs. The installed Python package separately supplies `pyalps::runtime` for extensions that share ALPS objects or HDF5 handles with pyalps. A matching C++ SDK is still required for headers and compile settings.

```cmake
find_package(Python 3.11 REQUIRED COMPONENTS Interpreter Development.Module)
execute_process(
  COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
  OUTPUT_VARIABLE nanobind_ROOT
  OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ERROR_IS_FATAL ANY)
find_package(nanobind 2.15.0 EXACT CONFIG REQUIRED)
execute_process(
  COMMAND "${Python_EXECUTABLE}" -m pyalps --cmake-dir
  OUTPUT_VARIABLE pyalps_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ERROR_IS_FATAL ANY)
find_package(pyalps CONFIG REQUIRED)
nanobind_add_module(my_module NB_STATIC my_module.cpp)
target_link_libraries(my_module PRIVATE pyalps::runtime)
```

`pyalps.get_cmake_dir()` exposes the same directory to Python tools. The C++ SDK neither installs this package nor discovers Python. The former SDK function `alps_target_link_pyalps` has been removed.

The same target supplies `<pyalps/export_simulation.hpp>` for exporting a derived simulation through nanobind. This Python-owned header replaces the old SDK header `<alps/ngs/detail/export_sim_to_python.hpp>`; update that include when rebuilding a downstream extension.

Wheel installation writes `pyalps/runtime.json`. After auditwheel or delocate repair, regenerate it with `python python/pyalps/_build_support/runtime_manifest.py --wheel path/to/pyalps.whl`. Cibuildwheel runs this automatically. The manifest records the final relative library paths; pyalps exposes these as imported CMake targets. On macOS, wheel finalization sets linkable `@rpath` library IDs and refreshes their signatures, so downstream builds need no binary-patching commands.

CMake derives build-time search paths from the imported targets. If you install or redistribute your extension, set its `INSTALL_RPATH` for the destination layout using normal CMake installation rules. The target does not hard-code the build environment's Python installation into installed extensions. For an extension installed for the same environment, CMake's `INSTALL_RPATH_USE_LINK_PATH` target property can retain the runtime search paths.

## Native Monte Carlo measurements

Python subclasses of `ngs.mcbase` use the same ALEA batches as native C++
simulations. Assign an accumulator to each named measurement, then append samples:

```python
from pyalps import alea, ngs

# Inside a simulation constructor, after super().__init__(parameters):
self.measurements["Energy"] = alea.BatchAccumulator(num_batches=64)
self.measurements["Correlations"] = alea.BatchAccumulator(size=length, num_batches=64)

# Inside measure():
self.measurements["Energy"] << energy
self.measurements["Correlations"] << correlations
```

`ngs.collectResults(simulation)` returns an owning `dict` of `BatchResult`
snapshots. Means and errors are vectors, including one-component scalar results.
Retrieved measurement handles remain usable after collection entries are removed
or replaced. `ngs.saveResults(results, simulation.parameters, archive, path)`
stores native kind-5 analysis results; `mcbase.save` stores complete kind-6 batch
state and the random stream. Subclasses save their own configuration and progress
alongside that base state. See the [Python simulation tutorial](../../tutorials/10-ngs/7_python_extend/).
The former NGS observable/result facades and experimental accumulator module are
removed. Released bins can be converted to native analysis results with the
[offline converter](../../src/tools/hdf5/README.md); analysis results do not supply
the sampling state needed for restart.

## Observable migration

The legacy `RealObservable`, `RealVectorObservable`, `RealTimeSeriesObservable`
and `RealVectorTimeSeriesObservable` producers are removed. Use native
`BatchAccumulator(size, num_batches)` for bounded weighted histories and
`AutocorrelationAccumulator(size)` for error-versus-batch-size diagnostics.
Use `size=1` for scalars. Snapshot an accumulator with `result()` before reading
its `mean`, `error`, `variance` or `tau`, and save through
`result.save(archive, path)`; the archive path now supplies the observable name.
`BatchResult.variance` describes batch observations, whereas
`AutocorrelationResult.level(0).variance` describes the original observations
when the accumulator uses its default `batch_size=1`.

`AutocorrelationResult.converged_errors` is an integer vector: 0 means a plateau,
1 means undetermined, and 2 means rising errors. It compares the last four
levels with at least 1024 batches each. Earlier errors below 90% of the final
error give status 1; below 82.4% give status 2. The strongest failed comparison
wins, correcting the old heuristic's ability to erase a failure when a later
comparison passes. Nonfinite errors or insufficient levels are undetermined.
This heuristic does not establish equilibration or ergodicity; constant data can
have a flat zero-error curve without exploring the state space. The complete
level hierarchy remains available, and `tau_available` reports sample sufficiency
rather than convergence.

Full chronological history uses ordinary NumPy arrays and HDF5 datasets, avoiding
a second observable implementation:

```python
import numpy as np
from pyalps import alea, hdf5

samples = np.asarray(samples, dtype=np.float64)  # (time,) or (time, component)
accumulator = alea.BatchAccumulator(1 if samples.ndim == 1 else samples.shape[1])
for sample in samples:
    accumulator << sample
with hdf5.archive("measurements.h5", "w") as archive:
    archive["/samples"] = samples
    accumulator.result().save(archive, "/simulation/results/Energy")
```

Accumulate those same samples into the native estimator while generating them.
For histories too large for memory, use an extendible h5py dataset. Never treat
averaged/compressed bins as individual chronological observations. Released
legacy observable files remain inputs to the offline converter; new producers
write native results and no longer emit legacy statistical schemas.

## MCData removal

`MCScalarData`, `MCVectorData` and the `pymcdata_c` extension are removed.
Measurement loading requires explicit conversion of legacy statistical records.
`--alea-results /simulation/results` converts a whole task file and prints the
choice for each observable: `--alea-batches` for recoverable linear/jackknife
histories, or `--alea-summary` to retain published statistics without
inferring history.
The primitive-normalization-only `--alea` option does not produce a modern
analysis result. Unconverted statistics now raise an actionable error rather
than being loaded implicitly or skipped after a logged exception.

Native results retain batches, covariance and diagnostics. Use their `transform`
method for correlated or nonlinear propagation. `FloatWithError` remains the
independent-error arithmetic tool for scalar or NumPy-vector summaries; Python 3
true division and correct unary negation are covered by the migrated arithmetic
checks. Plain deterministic eigenstate measurements and histograms remain
readable without statistical conversion. Lattice shape is no longer guessed
from `L` when reading a statistical vector; reshape components explicitly when
that matches the recorded geometry.

`loadBinningAnalysis` also requires native diagnostics. Convert released log-binning
hierarchies with `--alea-autocorr`; the converter recovers partial-bin weights and
preserves the source moments. The plotting reader now exposes every native level
and closes files on failure instead of logging an error and returning partial data.

## Merging and saving measurements

`DataSet.from_result(result, x=None, props=None)` creates plotting values while
retaining the native ALEA estimate in `native_result`. `loadMeasurements`
uses the same projection. Complex elliptic uncertainties are shown as circular
magnitudes; the full covariance remains available on the native result.

`mergeDataSets` pools independent runs of the same estimator family and component
coordinates, without mutating its input list or results. Use it only for runs
of the same physical ensemble. `mergeMeasurements` groups by observable name;
it does not decide whether runs with different parameters may be pooled.
`mergeMeasurementsFromFiles` now defaults to the native `/simulation/results`
path; pass `respath` explicitly for another layout.
Batches retain their sums and counts; summary moments retain their weights and
covariance. Autocorrelation merging pools the levels available in every nonempty
run, without treating independent runs as one continuous trajectory.

`saveMeasurements` writes those native results and component labels, encoding
observable names and preserving unrelated file contents. It does not infer
parameters from plotting properties; save parameters explicitly. Plot-only
summaries cannot supply the missing statistical evidence. Legacy files must
first be converted offline using the appropriate `alps-hdf5-convert` profile.

If you edit `DataSet.y`, merge/save will reject a mismatch with `native_result`.
For derived quantities, transform the result itself and construct a new dataset:

```python
import pyalps
from pyalps.dataset import DataSet

# For a joint BatchResult, propagate correlations through the jackknife.
ratio = joint_result.transform(lambda x: x[:1] / x[1:2])
data = DataSet.from_result(ratio, props={"observable": "Ratio"})
pyalps.saveMeasurements([data], "ratios.h5")
```

`alea.ReportedEstimate` is an explicit record for published statistics that lack
recoverable sampling evidence. Create it with `count=...`, a component-vector
`mean`, and optional `error`, `variance`, `tau` and `converged_errors` vectors.
Absent statistics remain absent attributes. It supports `read`/`save` and the
same `DataSet.from_result` plotting projection, but has no accumulator, covariance,
rebinning or native-pooling operations. `loadMeasurements` retains it in the
existing `native_result` field; inspect its type before requesting native-only
operations. Convert old files with `--alea-summary`; use `--alea-batches` instead
when complete histories can support native statistical analysis. Archived source
history remains in the converted file and is not copied into an in-memory report.

## Time-series analysis migration

`alea.mean`, `variance`, `error`, `autocorrelation`, the running means, cuts and
`make_dataset` now take real NumPy arrays with shape `(time,)` or
`(time, component)`. The `MCScalarTimeseries`, `MCVectorTimeseries` and view
classes are removed; cuts return NumPy slices. They accept read-only and strided
arrays. `ValueWithError` is also removed; use native results for statistical
estimates, or `pyalps.FloatWithError` for independent-error scalar arithmetic.

Autocorrelation returns an array of positive lags, starting at lag 1. It retains
sample-variance normalization: the denominator at lag k is `(N-k)*variance`,
with `variance` using `ddof=1`. Constant components have undefined correlation
and raise `ValueError`. Exponential fits return an ordinary `(amplitude, exponent)`
tuple instead of `StdPairDouble`. Fits use correct least squares in log space,
fixing the old regression's off-by-one sample count. Inclusive `from/to` bounds
refer to absolute one-based lags; negative bounds add the correlation length.
Threshold fits exclude the lower crossing. Invalid requests now raise exceptions
instead of printing usage or silently discarding nonpositive fit samples.

`integrated_autocorrelation_time` retains the positive-lag sum convention,
without adding 1/2. With a fit, it adds the continuum tail integral starting at
`N+1/2`; the variance inflation factor is `1+2*tau`. Binning errors use the native
ALEA level-selection policy. Inspect `AutocorrelationAccumulator` results for
level statistics and convergence; a returned error alone is not proof of convergence.

The C++ and Python examples in `tutorials/00-examples/alea` use
`generate_samples.py` to produce `timeseries.h5` with chronological AR(1)
observations. The former compressed-bin fixture has been removed. Offline
conversion preserves stored information but cannot reconstruct observations
lost to binning.

## HDF5 IO

`pyalps.hdf5.archive` owns an h5py file. Primitive datasets and attributes follow
h5py's dtype and shape rules; use NumPy arrays with explicit dtypes for scientific
fields. Reading a group returns an h5py Group. Dictionaries and ragged lists are
not inferred as containers: write their named fields explicitly. Modes are exactly
`r` (read), `a` (create/update), and `w` (truncate).

Scientific objects such as params, batch accumulators/results, RNGs, and simulations
retain their native save/load methods. These methods transfer file ownership for
the complete operation, closing h5py and reopening it afterwards. Native virtual
callbacks receive a small native archive that accepts scalars and explicit NumPy
arrays. Previously borrowed h5py groups/datasets become invalid during this
transfer; retained native callback views are closed before h5py reopens. Native
and h5py may use different HDF5 libraries and never share raw
HDF5 identifiers. Downstream C++ bindings that accept `alps::hdf5::archive&` use
`with archive.native() as native_archive:` around the complete native operation.

## Typed params and TOML migration

`ngs.params` now owns scalar and one-dimensional homogeneous values. Assignment
copies Python/NumPy input; retrieval returns a detached value. Reassign an edited
array (`p["x"] = values`) or use augmented assignment (`p["x"] += 1`). Missing
keys raise `KeyError`. Numeric strings, mixed Boolean/numeric arrays, arbitrary
objects, nested dictionaries, multidimensional arrays and `None` are rejected.
Python integers use signed 64-bit storage; integer-to-real conversion rejects
loss of precision. Boolean flags must use `True`/`False`.

Params checkpoints are explicitly versioned as `alps.params.v2`; old checkpoints
are not accepted by `ngs.params.load`. Use the [offline converter](../../src/tools/hdf5/README.md)
with `--parameters GROUP` to upgrade official flat parameter groups. The
unreleased `alps.params.v1` format is unsupported. The
analysis loaders also read the flat parameter groups of released result files. Checkpoints reject names and
string values containing NUL before overwriting stored parameters.

Installed C++ SDK consumers require the Boost version used to build that SDK.
The SDK exports the matching runtime search path on macOS, where Boost library
names alone may not distinguish incompatible versions from different providers.

MaxEnt, CT-HYB and CT-INT expose a common configured interface. MaxEnt's CLI
also accepts a TOML run file; see the [MaxEnt guide](../../src/apps/maxent/README.md).
The previous combined parameter dictionary, `AnalyticContinuation` function and
HDF5-as-run-file CLI are removed.

```python
from pyalps import cthyb, run_config
from pyalps.run_io import write_run_file

run = run_config.load("run.toml", cthyb.schema())
cthyb.solve(run)
write_run_file("copy.toml", cthyb.schema(), parameters=dict(run.parameters),
               input=dict(run.input), output=dict(run.output), execution=dict(run.execution))
```

For programmatic runs, `solver.prepare(parameters, input={}, output={},
execution={})` returns a validated `RunConfiguration` for `solver.solve(run)`. Native schema rules, defaults and application
checks apply to both file and programmatic runs. `output.results` is explicit.
Programmatic solver preparation anchors paths to the current working directory;
prepared paths remain stable if that directory changes. Input and output paths
in files resolve relative to the run file. Scientific
HDF5 datasets remain usable as input; their containers do not become run files.

`write_run_files(prefix, runs, schema=None, baseseed=None)` writes one TOML file
per explicit run plus a job manifest listing them; with a schema, the runs are
validated and the manifest names the application. Arrays remain values within a
run. The writer keeps supplied dictionaries unchanged, inserts generated seeds
in `execution.seed`, and refuses overwrites and colliding result/checkpoint
paths. When a schema is supplied, it also checks all declared input/output paths
across the job. Each active text-output directory is reserved for its run; other
runs' input and output paths must stay outside it.
`pyalps.run_io.execute(application, runs, mpi=None, mpirun="mpirun", concurrency=1)`
obtains each run's schema, rejects outputs that replace any run's inputs or other
outputs, and validates every run file or manifest with the application executable.
It then starts the runs in order, at most `concurrency` at a time, and returns their
absolute result paths in input order. After a failure it starts no further run and
raises once the active runs finish. With `mpi`, each run gets its own launcher;
pass `mpirun=["mpirun", "--bind-to", "none"]` so that concurrent MPI runs do not
bind to the same cores. `simplemc` now uses this native TOML workflow for all three
classical spin models, independent chains, exact HDF5 restart and direct VTK
snapshots; see the [simplemc guide](../../src/apps/mc/simple/README.md).
`spinmc` uses the same driver for Ising, XY, Heisenberg, O(4) and Potts models,
typed matrix arrays, local/cluster updates and native joint jackknife estimates;
see the [spinmc guide](../../src/apps/mc/spins/README.md). Its HDF5 results include
the derived statistics directly, so the separate spin evaluator is removed.
See the [CT-HYB guide](../../src/apps/dmft/qmc/hybridization/README.md)
for numerical formats and supported measurements.

## Installed-wheel validation

From the repository root, validate the installed package and retain reproducible
source, dependency, extension and distribution evidence:

```sh
python .github/scripts/validate_pyalps.py --output _build/validation --wheelhouse wheelhouse
```

The output includes logs, JUnit results, hashes, commands, timings and pass/fail
status. Failed runs retain their evidence. `--packaging` adds the release and
packaging tests. `--downstream` enables compiled SDK/extension consumers and
requires the matching SDK, CMake, a C++ compiler and nanobind. Check reported
skips; the ordinary run does not enable those consumers or provide mpi4py.

`--applications` additionally exercises input generation, dispatch and result
loading through six installed solvers. The four-site Heisenberg diagonalization
checks require ground-state energy -2 within 1e-10; short Monte Carlo and DMRG
runs check finite results without asserting convergence. Outputs stay beneath
the evidence directory. These workflows use the existing SDK and Python
environment without another build matrix.

MPI checks can run against the installed package separately:

```sh
mpiexec -n 2 python -m pytest -q \
  tests/pyalps/test_binding_surface.py::test_mpi4py_compatibility_surface \
  tests/pyalps/test_mpi_requests.py
```
