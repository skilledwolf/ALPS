# Changelog

User-facing changes and migration notes are recorded here, starting with the build modernization. Earlier releases are not yet catalogued in this file.

## Unreleased

### Changed

- Add `alps-hdf5-convert --core-parameters GROUP` for ALPSCore 2.3.3 `params` and
  dictionary checkpoints. It restores Core's escaped key names, reads its signed-byte
  Booleans and drops INI text, origins, help and descriptions, which belong to TOML
  run schemas. A fixture from compiled Core writers covers the profile, and the
  reconciliation probe's new Linux baseline checks that converted Core checkpoints
  load in ALPS. `src/alps/README.md` maps the ALPSCore params, archive and ALEA APIs
  to ALPS.
- Rename the older `alps::Parameters` to `alps::SymbolTable`
  (`<alps/expression/symbol_table.h>`) and reduce it to the in-memory text symbols of
  the expression, lattice and model libraries, owned by the `expression` module. Its
  text grammar, environment-variable substitution, XML, Osiris dump, Boost
  serialization and flat HDF5 input and output are removed, with `alps::Parameter`,
  `alps::ParameterList`, their parsers and the `legacy_parameters` module and its
  include cycle with `expression`. `alps::plot::Plot` takes `alps::params`.
- Read TOML run files in `checksign`, `printgraph` and `lattice2xml`, the model
  examples and the `08-alpsize` parameter and lattice lessons. The lattice and model
  library tests set their parameters in code. Nothing reads the older parameter
  text grammar any longer apart from its own tests.
- Convert typed parameters to the lattice and model libraries in one place.
  `graph_helper` and `model_helper` accept `alps::params`, `alps::seed_disorder`
  replaces five application copies of the disorder seeding, and `ALPS::params` no
  longer depends on the older parameter classes.
- Port classical energy `wanglandau` learning, overlapping-window weights,
  microcanonical sampling and temperature reweighting to TOML and native ALEA.
  Preserve histogram/round-trip diagnostics, native continuation and MPI chains;
  use joint batches for correlated thermodynamic errors. Include both energy
  endpoints in flatness checks, distinguish inaccessible bins from zero log
  weights, and require reference normalization for absolute free energy/entropy.
  Remove the replaced Parapack Wang–Landau/XDR implementation and fixtures.

- Port the `ising` and `heisenberg` library examples to TOML and native ALEA,
  preserving bond-type couplings, normalized magnetization moments, both
  Heisenberg Binder ratios and replica exchange. Share ladder schemas and
  feedback with `loop`; retain the classical infinite-temperature endpoint.
  Correct self-bond update weights: those bonds contribute constant energy.
  Native checkpoints support partial-run continuation, production extension
  and MPI redistribution of independent ladders. Remove the replaced workers
  and their obsolete scheduler output fixtures.

- Distribute the walkers of one replica ladder over MPI ranks with
  `execution.parallel = "replicas"` in `loop` and the `ising` and `heisenberg`
  examples; the default `"chains"` distributes independent ladders. Each
  temperature keeps one chronological measurement stream, and checkpoints
  resume with either layout and any process count.

- Port the `exchange` library example to TOML and native ALEA. It accepts the
  former `"ising; exchange"` and `"multiple parallel ising; exchange"`
  algorithms; `execution.processes_per_walker` (formerly `PROCESS_PER_WORKER`)
  decomposes each ring walker over a team of MPI processes. Every layout
  reproduces the serial ladder exactly. Run `"loop; exchange"` jobs with the
  `loop` application, which also replaces the `loop_single` example program.

- Remove the Parapack framework (`<alps/parapack/*.h>`) and its `pevaluate` and
  `poutput` commands; every former Parapack example now runs on the native
  Monte Carlo driver and moves from `tutorials/00-examples/parapack` to
  `tutorials/00-examples/mc`, keeping its program name. The SDK no longer links OpenMP: `ALPS_ENABLE_OPENMP`
  only builds the OpenMP thread-safety tests, and the classical examples take
  their own `ALPS_ENABLE_OPENMP` option. Parapack's OpenMP threads per clone have
  no replacement; chains run as MPI processes. See [migration](#removed-and-migration).

- Remove the scheduler's Monte Carlo layer (`<alps/scheduler/montecarlo.h>`)
  and the `compactrun` command for its run dumps. `convert2xml` keeps its
  parameter-file and spectrum conversions and rejects Monte Carlo XDR
  checkpoints and XML runs with the migration path instead of converting them.
  The DMRG and diagonalization tasks keep the generic scheduler.

- Remove the legacy ALEA observables (`<alps/alea.h>`, `ObservableSet`, the
  binning strategies and evaluators) from `ALPS::alps` after their last
  consumers migrated. `ALPS::statistics` is the only statistics library; the
  `fleas_correlated` sampling example now uses its autocorrelation estimator.

- Store native autocorrelation levels as one dataset per field, levels first
  (`level/count`, `level/count2`, `level/mean`, `level/var` in results and the
  matching `levels/` fields in checkpoints) instead of one group per level.
  An observable now takes about 9 HDF5 objects per chain instead of about 60,
  which halves the time a short simulation spends writing its outputs.
  `alps-hdf5-convert --core-alea autocorr` and `--alea-autocorr` write the
  stacked levels.

- Remove the unbuilt `src/apps/qmc/sse` and `sse2` sources, which the released
  QMC build already excluded and which no longer compiled. The installed
  `dirloop_sse` is unchanged. The `sse2` optimized-ensemble workflow
  (`SIMULATION_PHASE=2` with external `LOGG_FILENAME` weights) has no
  replacement; `qwl` provides Wang–Landau and multicanonical sampling.

- Run several TOML runs at once with
  `pyalps.run_io.execute(application, runs, concurrency=N)`, restoring the
  former schedulers' concurrent execution of a job's tasks. Every run is still
  validated before the first starts, and no further run starts after a failure.
  `mpirun` also accepts a launcher argument list, such as
  `["mpirun", "--bind-to", "none"]` for concurrent MPI runs.

- Port the Fortran bridge (`ALPS::fortran`) to TOML runs, native ALEA and typed
  HDF5 checkpoints. Fortran callbacks now use `bind(C)` interfaces with an explicit
  caller handle; see [migration](#removed-and-migration).

- Port `loop` to native ALEA and the TOML CLI, retaining continuous-time/SSE,
  annealing, signed/improved/custom estimators, replica exchange and both
  temperature-feedback methods. Checkpoint physical walkers, exchange state,
  feedback and unfinished batches; keep results separate by temperature.
  Correct normal-estimator sign tracking after cluster flips and include
  longitudinal-field energy in replica exchange weights. Retain one common
  Hamiltonian across the ladder, including parameter-dependent XML couplings.
  Finish released simulations with ALPS 3.0; the converter migrates results for analysis.

- Port `worm` and `dirloop_sse` to native ALEA and the established TOML CLI,
  retaining executable names, update kernels and both RNGs. Checkpoint complete
  worldline/operator state, tuning and unfinished batches; retain per-chain
  evidence and sign-aware joint analysis. Correct Green-function signs and
  interval normalization, custom correlation normalization, chain-density
  moments and centered compressibility. Fixed-number worm tuning preserves
  sector weights and reports energies for the requested Hamiltonian. Reject
  unimplemented estimator combinations. Migrate the Python tutorials and English/Japanese notebooks to TOML and native
  results; remove obsolete XML input examples. Finish released runs with ALPS 3.0; see the quantum solver run guide.

- Port `qwl` to native ALEA and typed TOML runs using the shared Monte Carlo
  runner. Preserve refinement modes, magnetic measurements, intermediate
  coefficients, serial/MPI chains and both RNGs; add exact operator/histogram
  restart and production extension. Correct expansion-window initialization,
  combinatorial coefficient normalization and independent-chain evaluation.
  Migrate the evaluator and tutorials; provide offline conversion of released
  per-run final estimates with `alps-hdf5-convert --qwl-sites`.

- Port `spinmc` to typed TOML runs and native ALEA for Ising, XY, Heisenberg,
  O(4) and Potts models, preserving local and legal cluster updates, matrix
  couplings, onsite interactions and fields. Retain exact spin/RNG/progress and
  partial-bin checkpoints; pool independent raw chains before joint nonlinear
  analysis. Keep direct moments separate from improved cluster estimators.
  Correct field-aware cluster selection, signed/mixed matrix handling,
  asymmetric-bond and self-loop local deltas, onsite constants and staggered
  cluster projections. Remove the old scheduler/factory/matrix kernels and
  separate evaluator; migrate classical tutorials and notebooks. Chains run
  in serial or MPI with `mt19937` or `lagged_fibonacci607`; legacy XML runs reject.
- Share private native Monte Carlo run orchestration and vector-aware raw-chain
  pooling across applications. Keep exact restart and output validation in one
  path. Omit optional thermodynamic estimates when moment subtraction cannot
  resolve their fluctuations, including every jackknife estimate; preserve raw
  evidence instead of amplifying roundoff at very low temperatures.
- Center native ALEA batch sums before calculating variance and covariance.
  Constant measurements and small fluctuations around large means now retain
  finite uncertainties without subtracting nearly equal raw second moments.
  Preserve unequal and partial batch weights and complex covariance conventions.
- Port `simplemc` to typed TOML runs and native ALEA batches, sharing one local
  update engine for Ising, XY and Heisenberg models. Preserve custom lattices,
  graphs, bond couplings and magnetic fields; exclude thermalization and pool
  independent raw chains before nonlinear analysis. Retain exact RNG, spin and
  partial-bin restart state, validate all tasks before execution, and emit VTK
  snapshots directly. Correct self-loop update weights and zero-warmup sampling.
  Remove the duplicate Parapack workers and evaluators; update the snapshot
  tutorial. Independent chains run in serial or MPI with `mt19937` or
  `lagged_fibonacci607`; legacy XML inputs are rejected.
- Validate native ALEA checkpoints by replaying the shared batch merge algorithm.
  Reject inconsistent cursors, counts, offsets and nonzero empty bins before
  replacement; remove manual modular cursor reconstruction and duplicate
  application validation. Keep batch sampling explicitly unit-weighted.
- Correct native ALEA Hotelling mean tests: retain full covariance and fractional
  effective counts, pool two-sample degrees of freedom correctly, and test complex
  batches as joint real/imaginary data. Reject missing covariance and insufficient
  observations. Use Boost.Math F tails and standard upper-tail p-values; remove
  synthetic pooled-result state and the custom special-function implementation.
- Remove unbuilt spin-fitting code, unused TCP/ULFM scheduler remnants and dead
  checkpoint branches. Correct scalar summary pooling and empty-simulation
  measurement collection while the remaining scheduler applications migrate.
- Move NGS `mcbase`, its MPI adapter and Python simulation bindings to native ALEA
  batches. Use explicit component dimensions and owning result snapshots; retain
  complete RNG and unfinished-bin checkpoint state, with staged base loads.
  Include empty ranks in collective result reduction and validate matching result
  requests before pooling. Propagate scheduled rank-local failures. Remove the
  old measurement/result facade, scalar feature-stack API, five Python statistical
  modules and unused scheduler prototypes. Consolidate duplicate Ising tutorials,
  repair virtual checkpoint hooks, correct Heisenberg coupling acceptance and
  remove its invalid partial-magnetization susceptibility estimator.
- Add `--alea-batches GROUP` to the offline converter. Recover released ALPS 3.0.0
  complete linear bin sums and partial weights as native ALEA analysis results;
  recompute uncertainty and preserve vector covariance. Reject unrecoverable
  histories instead of inventing restart state or missing covariance.
  `--alea-results GROUP` applies this to every observable of a released task
  file and keeps reported estimates for observables without bins. The total-sum
  check uses the writer's summation round-off bound, so genuine released files
  with accumulated rounding convert.

- Apply opposite magnetic-field signs to paired flavors in the general
  paramagnetic DMFT Hilbert transform; previously every flavor used `MU - H`.
- Port general multiband CT-INT density interactions to `interaction`, using one
  kernel for two-flavor and multiband runs. Correct sparse-matrix detailed balance
  with uniform unordered-pair proposals; isolated flavors and zero interactions
  terminate normally. Validate finite, symmetric, zero-diagonal matrices, expand
  flavor-dependent TOML/Python schemas, and retain full native ALEA density moments
  for general Fourier tails. DMFT now uses external solver processes exclusively;
  remove its `Interaction Expansion` selector, scheduler adapters and duplicate
  CT-INT implementation. Update tutorials to select `interaction` and remove
  `output.matrix_size`; order diagnostics remain in canonical HDF5 results.
- Migrate Hirsch-Fye to native ALEA joint Green/sign batches and raw MPI
  collection. Form physical endpoints before accumulation, retain partial bins,
  exclude warm-up, and publish native statistics, Green functions and provenance
  atomically. DMFT now selects the external `hirschfye` executable; remove the
  `Hirsch-Fye` scheduler selector, nonfunctional restart methods and disabled
  four-point code/settings. Preserve its seeded generator sequence and update
  cadence, check complete matrices and LAPACK solve failures, bound the rebuild
  interval and correct the initial determinant sign and negative-weight
  heat-bath acceptance. Propagate preparation, sampling and output failures
  across MPI ranks.
- Migrate CT-HYB to native ALEA batches and MPI collection. Normalize physical
  measurements with the sign from their own measurement cadence, pool raw
  replicas before analysis, and retain G2/H2 component errors in linear memory.
  Publish native kind-5/kind-2 results, derived Green functions and provenance
  together; Python analysis reads both families. Form time endpoints before
  accumulation so means, errors and covariance agree, correct F covariance
  selection and the missing sign in the zero-frequency density correlator.
  Replace `NUM_BINS` with even `execution.bins` and remove `ACCURATE_COVARIANCE`.
  Result files remain analysis outputs without a solver restart interface.
- Add componentwise real/imaginary ratio propagation for native elliptic ALEA
  results, retaining numerator/denominator covariance with linear storage.
  Expose squared weights and ordinary standard errors through the existing
  heterogeneous result facade; reject invalid ratio domains without mutation.
- Migrate standalone CT-INT to native ALEA batches and MPI reduction. Normalize
  sign-weighted W and density estimates using joint numerator/sign batches;
  retain partial replica bins and unsigned diagnostics. Write versioned batch
  results and publish them with Green functions and run provenance atomically.
  Python measurement loading retains modern result errors. Local batch slots
  must be even and at least two; these analysis outputs are not solver checkpoints.
  Preserve the chain sign during inverse-matrix rebuilding and remove the
  implicit, incorrect `staggered_sz` text output; measured spin results remain.
  Derive zero Fourier moments for the analytic atomic input and reject
  contradictory explicit moments, correcting the noninteracting histogram output.
  Use fermion occupation idempotence for same-site density products and spin
  squares instead of squaring conditional density estimates.
- Keep empty runs neutral during modern autocorrelation reduction, preserving
  the depth and errors of contributing runs while rejecting malformed levels.

- Preserve independent-run ALEA bins and partial-bin weights during reduction;
  retain only autocorrelation levels represented in every run. Stage reductions
  before publishing results and reject malformed shapes and weights collectively.
  Correct elliptic complex variance centering.
  Add an optional native MPI reducer with unsigned 64-bit sample counts; custom
  reducers must implement the unsigned overload and downstream binaries need a
  rebuild. Remove the unused unsigned 32-bit reinterpretation overload.
- Correct signed ratio uncertainty propagation. Use the existing joint batches
  and transforms, with full covariance,
  relative central differences and weighted jackknife; reject singular domains
  and malformed bins, and preserve tiny-bin contributions for linear transforms.
- Remove unused Boost tuple/shared-array HDF5 adapters and their test machinery;
  the archive-copy test uses a standard owning buffer.
- Expose modern real/complex batch ALEA in Python and migrate the pure Python
  Ising tutorial to the shared native checkpoint codec and publication helper.
  Reject malformed RNG checkpoints without changing the stream. Fix squared
  batch-weight overflow and reject wrong-sized samples before changing state;
  remove the unsafe autocorrelation accumulator/result merge.
- Remove the unsupported `ALPS_NGS_USE_NEW_ALEA` backend selector, unused wrappers
  and unregistered tests; the remaining scalar feature-stack API is now retired.

- Consolidate ALPSCore's Eigen-based modern ALEA as independent `ALPS::statistics`,
  with one canonical HDF5 adapter and explicit result versions/kinds. Add actual
  resumable batch-accumulator checkpoints and migrate the accumulator-only Ising
  tutorial, preserving RNG and unfinished batch state. Fix incoming Eigen layout
  serialization and nonsquare linear-transform dimensions.
- Replace legacy statistical objects only after successful load; remove stale
  optional datasets on save. Publish scheduler HDF5 snapshots through checked
  close while preserving their existing XDR/XML coordination.
- Remove dormant deprecated accumulator code and its obsolete include exemptions.
  Convert released ALPSCore 2.3.3 ALEA results offline, with a compiled Core fixture
  and native readers covering covariance and per-batch counts.

- Replace the hand-written HDF5 backend with HighFive 3.3.0 and require HDF5 2.x. Complex numbers use ordinary compound datatypes, Booleans use enums, and empty arrays keep their zero-length dimensions. HighFive is private to the implementation; installed SDK consumers do not need it. Remove legacy complex markers, numeric/string coercions and old pair/matrix readers from the runtime.
- Add the standalone `alps-hdf5-convert` tool (Python, h5py and NumPy) to translate legacy complex and marked Boolean encodings into ordinary HDF5 compound/enum datatypes. Explicit `--parameters GROUP` and `--alea GROUP` profiles migrate the flat parameter and observable/result schemas written by ALPS 3.0.0. It writes a separate file and rejects ambiguous type/shape reconstruction. See the [conversion guide](src/tools/hdf5/README.md).
- Simplify archive ownership: independent opens have independent permissions, modes are exactly `r` (read), `a` (create/update) and `w` (truncate), and explicit close invalidates copied views. Remove family, memory and compression mode modifiers and the global filename registry. Publish full checkpoints through explicit `alps::hdf5::save_checkpoint`, preserving the previous file on serialization failure and closing the temporary file before publication.
- Use h5py for ordinary Python HDF5 data, attributes and groups. NumPy/h5py determine dtype and shape; Python lists and groups no longer use ALPS-specific container inference. Native scientific serializers transfer file ownership for a complete operation through `archive.native()`. Python virtual checkpoint callbacks receive an owned native view with scalar/NumPy field IO; views become closed when the operation ends.
- Consolidate typed params around ALPSCore's owning dictionary/value model (reference `7146b9e1f017938a94e5dae35d88467cc5ba7969`). Store fixed-width integer families, real/complex scalars and homogeneous vectors; reject lossy or string-to-number conversions. TOML loading, schema defaults/validation and run provenance belong to the new `ALPS::run_config` component, using toml++ ≥ 3.4 privately.
- Migrate MaxEnt to TOML run files and a separate application schema, scientific-data input and explicit output/execution settings. C++ and Python use the same native validation. Preserve numerical algorithms and reference cases; store frequency-bin widths in the result instead of unconditionally writing `deltaOmega.dat`. See the [input guide](src/apps/maxent/README.md).
- Migrate CT-HYB (`hybridization`), CT-INT (`interaction`), Hirsch-Fye (`hirschfye`) and the DMFT driver (`dmft`) to TOML run files. Each application embeds its schema, prints it with `--schema` and installs it under `share/alps/schemas`; `--validate` checks a run without sampling or writing results. Run loading rejects path outputs that would replace the run file, an input or another output. The DMFT driver composes its schema with the selected solver's, launches external solvers from the ALPS `bin` directory and passes per-flavor `/Delta_<f>` or `/G0_<f>` vectors. See the [DMFT guide](src/apps/dmft/qmc/README.md) and the [CT-HYB guide](src/apps/dmft/qmc/hybridization/README.md).
- Migrate the exact diagonalization programs `sparsediag` and `fulldiag` to TOML run files. Both embed their schema, print it with `--schema`, install it under `share/alps/schemas`, check runs with `--validate` without diagonalizing, and accept several run files per invocation. Undeclared model parameters are accepted with inferred types; lattice and model libraries are `input.lattice_library` and `input.model_library`. Results keep the `/spectrum/sectors` layout read by `pyalps.loadSpectra` and `pyalps.loadEigenstateMeasurements`, with typed `/parameters`, `/run_config` and `/spectrum/number_of_sites`. `fulldiag_evaluate` reads these results without rebuilding the lattice and writes the same XML plots to `<prefix>.plot.*.xml`; `--h`/`--H` now evaluates the spectrum at that field, and `T_MIN` defaults to `DELTA_T`.
- Migrate `dmrg` to TOML run files with the same command line, validation and library inputs as the exact diagonalization programs, which completes the retirement of the XML scheduler applications. `STATES` is an integer array, and DMTK's temporary files go to `execution.temporary_directory` (by default the system temporary directory) and are removed when the run ends. Results keep the `/spectrum` and iteration-history layout read by `pyalps.loadEigenstateMeasurements` and `pyalps.loadMeasurements`, with typed `/parameters` and `/run_config`. The DMRG tutorials and English notebooks write run files and run under Python 3.
- Change params checkpoints to `alps.params.v2`, with native datatype/rank determining value types, native complex/Boolean datatypes, ranked empty vectors and transactional loading. Indexed name/value entries preserve arbitrary parameter names without duplicate type tags. Python params copy values on assignment and retrieval; missing keys raise `KeyError`. Native spin tutorials use TOML input with a shared application schema.

- Export `ALPS::configuration`, `ALPS::containers`, `ALPS::numerics`, `ALPS::numeric_io` and `ALPS::solver_headers` as interface targets. Foundation libraries and interfaces own their header sets without inheriting aggregate `ALPS::headers`. Separate array storage from mathematical helpers and numerical types from HDF5 adapters, removing the foundation include cycle; the expression/older-parameters cycle remains. Numerical algorithms and the moved headers' public include names are preserved, with no ALPSCore implementation imported. Params conversion headers move to `params/adapters/include/` and remain exposed through the aggregate interface with implementations in `ALPS::alps`.
- Group `src/tools/` commands by parameter preparation, lattice export, scheduler formats, diagnostics and XML transformations. Keep historical inactive tools under explicit owners without enabling them. Executable names, installation components and source contents are preserved; `pconfig` now links only the utilities component.
- Organize `src/alps/` by semantic ownership, including numerical/container helpers, XML, older parameters, models, observables and execution, replacing the transitional `common/` and `runtime/` groups. Modules own their headers, sources and local tests; configuration templates live in `cmake/config/`. Explicit file sets preserve public include names, and generated headers use `<build-dir>/generated/include/alps/`. A CMake-generated module manifest supports ownership and include-dependency checks. Remaining simulation modules use the aggregate header interface; physical owners are not all independent libraries.
- Export `ALPS::xml` for XML parsing/output and `ALPS::cli` for the existing `mcoptions` and `parseargs` command-line grammars. MaxEnt's executable now links `ALPS::cli` instead of `ALPS::alps`; its HDF5 input, command-line behavior and scientific calculations are unchanged. Move `<alps/plot.h>` to the `plotting` header module because it also uses older parameters. Public include names remain stable, and `ALPS::headers` remains the aggregate compile interface. Python runtime packaging includes both libraries; rebuild downstream binaries after the split.
- Export `ALPS::osiris` for dump/XDR serialization and process/communication state. MaxEnt's solver links the foundations, Osiris and numerical providers without the simulation runtime, preserving its numerical calculations and stop-callback behavior. Python runtime packaging includes Osiris. Rebuild downstream binaries after the library split.
- Export `ALPS::params` as a separate typed parameter library with component-owned exports. Isolate the existing parameter-file constructor, XML input and older `Parameters` conversion in adapters owned by `ALPS::alps`; their source APIs and parsing behavior are preserved. Python packages carry the shared runtime components. Rebuild downstream binaries after the split.
- Export `ALPS::hdf5` as a separate archive library with its own symbol exports. It links utilities and HDF5/Boost dependencies without the simulation runtime; `ALPS::alps` links it transitively. Python packages include one copy of this runtime component. NGS termination polling belongs to utilities; opening an archive no longer installs fatal-signal handlers. Rebuild downstream binaries after these ownership changes; public include paths are unchanged.
- Export `ALPS::utilities` as a separate library with its own symbol exports. `ALPS::alps` links it transitively, and Python packages carry this runtime component. Rebuild downstream binaries after this library split; source include paths and utility APIs are unchanged.
- Group utilities, HDF5 and NGS parameter headers, implementations and tests by module under `src/alps/`; separate MaxEnt's implementation, CLI and tests under `src/apps/maxent/`. Public include paths and exported library targets are unchanged. See the [module layout](src/alps/README.md) for the ALPSCore reconciliation boundaries.

- Require Python 3.11 or newer for pyalps, its wheels and the Python scripts of native tests.
- Require CMake 3.27 or newer and an externally installed Boost 1.76 or newer with CMake packages. The SDK requires C++17/C11 compilers, HDF5's C library, and LP64 BLAS/LAPACK; bundled Boost builds and alternate numerical integer/symbol ABIs are no longer supported.
- Export CMake targets for the installed SDK, applications and solver libraries. Downstream C++ projects link `ALPS::alps`; Python extensions sharing pyalps objects use `pyalps::runtime`. MaxEnt, CT-HYB and CT-INT Python wrappers link the SDK's solver libraries instead of compiling their implementations again.
- Make MPI opt-in with `ALPS_ENABLE_MPI=ON`. Standalone builds enable applications and native tests by default; embedded `add_subdirectory` builds default to the library alone. The default SDK uses shared libraries, as required by the Python bindings.
- Support editable Python development against an installed SDK, with packaged native runtime resources and explicit downstream ABI checks. Redistributable wheels bundle native dependencies through platform repair and are tested after installation on fresh runners.
- Organize tutorials into numbered topic directories, with introductory scripts directly in `tutorials/01-intro/`. `tutorials/00-examples/` is an independent API reference collection. Tutorial sources are installed only when the `tutorials` component is requested. See the [tutorial index](tutorials/README.md).

### Removed and migration

HDF5 archive construction uses exact string modes: replace integer `0`/`READ` with `"r"` and `1`/`WRITE` with `"a"`. Use `"w"` to create a fresh file, `"a"` to preserve other groups, and `alps::hdf5::save_checkpoint(filename, callback)` for atomic full checkpoint publication. Character and integer constructors, public property flags and mode modifiers are removed. Copies share a file owner; destruction releases a view, while explicit close ends the shared operation. Unsafe implicit copy assignment is deleted. The fatal-signal APIs `archive::abort()` and `signal::listen()`/`segfault()` are removed; the library leaves the application's crash handlers in place.

Archive `extent()` now reports actual HDF5 dimensions: a scalar returns `{}`, and a complex vector of length N returns `{N}`. Empty pointer writes require an explicit zero-length shape. `set_complex()` and the marker traits are removed. Numeric conversions use HDF5; string/numeric and Boolean/numeric coercions are rejected. Convert old files offline with `alps-hdf5-convert SOURCE DESTINATION` before reading them with the new runtime. The converter cannot infer lost shapes or undocumented application schemas; review its report before treating an old solver checkpoint as resumable.

Typed params no longer accept legacy text/XML files or old parameter checkpoints. Include `<alps/params.hpp>`, use `.exists()` / `.as<T>()` / `.value_or(key, fallback)`, and supply Boolean flags as actual Booleans. Native non-const `[]` inserts an unset entry; const lookup throws, and iteration is lexicographic. Convert official flat parameter groups with `alps-hdf5-convert SOURCE DESTINATION --parameters /parameters`; unmarked signed-byte Boolean intent still requires an explicit declaration. Unreleased v1 checkpoints are not a supported migration input.

Python archive indexing reads primitive datasets or returns h5py groups; use explicit group/dataset operations rather than assigning nested dictionaries or inferring lists from numbered children. Native object `save`/`load` operations close h5py temporarily, invalidating previously borrowed h5py objects. Deprecated `h5ar`, `iArchive`, `oArchive` and XML-export aliases are removed.

The unreferenced `alps::ngs_parapack` XML frontend is removed: all its references were internal to that frontend, and the project assumes no external consumers. The `alps::parapack` framework is also removed; write simulations against `alps::mc::main` as in the classical examples of `tutorials/00-examples`. Its `pevaluate` and `poutput` commands read released Parapack XDR and HDF5 results: convert results with `alps-hdf5-convert SOURCE DESTINATION --parameters /parameters --alea-results GROUP` and analyze them with pyalps, and finish released mid-run Parapack simulations with ALPS 3.0. The scheduler's Monte Carlo classes `MCRun`, `MCSimulation`, `LatticeMCRun`, `LatticeModelMCRun`, `DummyMCRun` and `SimpleMCFactory` are removed as well; `compactrun` and the Monte Carlo run and simulation conversions of `convert2xml` follow the same migration path. The legacy ALEA headers `<alps/alea.h>` and `<alps/alea/*.h>` are removed: accumulate with `ALPS::statistics` (`<alps/alea.hpp>`) and convert released observable files with the `alps-hdf5-convert` ALEA profiles. Older `alps::Parameters` applications remain; migrating their application orchestration is separate work. The lattice and model libraries take typed parameters directly: construct `alps::graph_helper` and `alps::model_helper` from `alps::params`, and call `alps::seed_disorder(p)` where an application seeded the disorder generator itself. `<alps/ngs/make_deprecated_parameters.hpp>` is removed; `explicit alps::SymbolTable(params const&)` replaces `make_deprecated_parameters`.

The scheduler library is removed with its last XML-job applications: `<alps/scheduler/...>` task, worker and dispatch classes, `DiagTask`, and the `convert2xml` and `parameter2xml` commands. Every application reads TOML run files and writes HDF5 results that pyalps loads directly. `snap2vtk` remains for snapshot files written by released versions.

`alps::mcmpiadapter` is removed with the communicator overload of `alps::stop_callback`: MPI simulations run on `alps::mc::main`, which spreads `execution.chains` independent chains over the ranks, stops them at `execution.time_limit` and pools their results. An `mcbase` simulation adds a `(parameters, bins, chain)` constructor, `checkpoint_parameters` and `completed_sweeps`, as in the [mcbase tutorial](tutorials/09-code/07-mcmain-mcbase/mpi.cpp). That tutorial's one MPI program also runs parameter scans and replaces the `10-ngs` MPI copies.

The Python MaxEnt, CT-HYB and CT-INT modules expose `schema()`, `prepare(parameters, input, output, execution)` and `solve(run)`. CT-INT's `schema(parameters={})` expands its per-flavor settings. `solve` takes a run returned by `prepare` or loaded with `pyalps.run_config.load` instead of a combined parameter dictionary, and MaxEnt's `AnalyticContinuation` function is removed. `pyalps.runDMFT` is replaced by `pyalps.run_io.execute(application, runs)`, which validates every TOML run file or job manifest before running any; `write_run_file` and `write_run_files` take an optional schema. DMFT's in-process solvers and scheduler settings `NRUNS`, `CONVERGENCE_CHECK_PERIOD` and `SWEEP_MULTIPLICATOR` are removed.

The old build interfaces are removed without compatibility aliases. Reconfigure existing build scripts and downstream projects using [CONTRIBUTING.md](CONTRIBUTING.md#getting-started-with-the-code).

Numerical persistence now requires an explicit adapter: include `<alps/hdf5/matrix.hpp>` or `<alps/hdf5/numeric_vector.hpp>` and link `ALPS::numeric_io`. The matrix umbrella `<alps/numeric/matrix.hpp>` no longer includes HDF5. HDF5 `save`/`load` members on diagonal matrices and deprecated BLAS matrix/vector classes, plus numerical matrix `write_xml` and XML insertion overloads, were unused in this repository and are removed without replacement adapters. Header ownership moves preserve public include spellings, but these method removals and explicit adapter requirements are source API changes.

| Previous interface or location | Replacement |
| --- | --- |
| `wheel-deps` preset | `distribution` preset; update `ALPS_DIR` to the new installation |
| `ALPS_BUILD_TESTS` | `ALPS_BUILD_TESTING`; a parent's `BUILD_TESTING` does not control ALPS tests |
| `ALPS_BUILD_LIBS_ONLY=ON` | `ALPS_BUILD_APPLICATIONS=OFF`; the `sdk` preset also disables tests |
| `ALPS_BUILD_EXAMPLES` | Build the standalone projects in `tutorials/00-examples/` against an installed SDK |
| `ALPS_BUILD_FORTRAN` | The SDK always includes `ALPS::fortran`; build Fortran examples separately with a Fortran compiler |
| `ALPS_INSTALL_HEADERS` | SDK headers are always installed |
| `ALPS_INCLUDE_TUTORIALS` | `cmake --install <build-dir> --component tutorials` |
| `ALPS_ENABLE_OPENMP_WORKER` | Removed with the Parapack worker; the classical examples in `tutorials/00-examples` take their own `ALPS_ENABLE_OPENMP` |
| `ALPS_USE_MKL_PARALLEL` | `BLA_VENDOR=Intel10_64lp` for threaded or `Intel10_64lp_seq` for sequential MKL |
| `ALPS_BUILD_ON_CRAY` | Removed without replacement; CMake's BLAS and LAPACK finders select the numerical libraries |
| `ALPS_LINK_BOOST_TEST`, `ALPS_INSTALL_BOOST_TEST` | Removed; native tests use header-only Boost.Test, which the SDK does not install |
| `ALPS_BUILD_DEVELOPER_TOOLS` | Removed without replacement |
| `ALPS_NGS_OPENMPI_ULFM` | Removed with the unused ULFM fault-tolerance sources |
| `PYALPS_BUILD_APPLICATIONS` | `PYALPS_BUILD_SOLVERS`; `PYALPS_BUNDLE_APPLICATIONS` still controls bundled executables separately |
| `ALPS_USE_SYSTEM_BOOST`, `Boost_ROOT_DIR` / bundled Boost discovery | External Boost packages are always required; locate them through `Boost_ROOT` or `CMAKE_PREFIX_PATH` |
| `LAPACK_64_BIT`, alternate `BIND_FORTRAN_*` ABIs | Use LP64 BLAS/LAPACK with lowercase, trailing-underscore symbols |
| `UseALPS.cmake`, `include.mk`, `alpsvars` scripts | Imported SDK targets and explicit installation paths; add the installed `bin` directory to `PATH` |
| `UsePyALPS.cmake`, `<alps/ngs/detail/export_sim_to_python.hpp>` | [pyalps downstream CMake package](python/pyalps/README.md#downstream-native-extensions) and `<pyalps/export_simulation.hpp>` |
| `plot2*`, `convert2html`, `convert2text`, `extract*` shell tools | `alps-xml plot`, `alps-xml convert`, `alps-xml extract`; see [XML tools](CONTRIBUTING.md#xml-resources-and-tools) for formats and dependencies |
| `xslttransform`, `use_local_stylesheet`, `archivecat` shell tools | `alps-xml` for rendering; run `xsltproc` with the stylesheets in `share/alps/xml` for other transformations |
| `archive` (SQLite result index), `xml2archive`, `txt2archive` | Removed without replacement; load results with pyalps |
| `spinmc_evaluate`, `loop --evaluate`, `pyalps.evaluateSpinMC`, `pyalps.evaluateLoop` | `spinmc` and `loop` write derived results at the end of each invocation; convert released results with `alps-hdf5-convert --alea-results` |
| `pyalps.runApplication` | `pyalps.run_io.execute(application, runs)`; every application reads TOML run files |
| `pyalps.writeInputFiles`, `writeTaskXMLFile`, `writeParameterFile`, `generateSeed`, `input2output` | `pyalps.run_io.write_run_files(prefix, runs)`; seeds are set in each run's `[execution]` section |
| `parameter2xml` | TOML run files, for example written with `pyalps.run_io.write_run_files(prefix, runs)` |
| `alps::Parameters` (`<alps/parameter.h>`) as the lattice and model symbols | `alps::SymbolTable` (`<alps/expression/symbol_table.h>`), built from `alps::params` |
| `alps::Parameters` text, XML, dump and HDF5 input/output; `alps::Parameter`, `alps::ParameterList` | `alps::params` with TOML run files and `alps.params.v2` checkpoints; convert released flat parameter groups with `alps-hdf5-convert --parameters GROUP` |
| `pyalps.writeInputH5Files`, `pyalps.getInputH5Files` | `pyalps.run_io.write_run_files(prefix, runs)`; `pyalps.getParameters` reads a results file's `/parameters` |
| `parameter2hdf5`, `p2h5` and their `.in.h5` parameter files | TOML run files; `pyalps.run_io.write_run_files(prefix, runs)` writes them |
| Parameter text on standard input or as a file for `checksign`, `printgraph` and `lattice2xml`; `checksign` XML jobs and `{ }` task lists | One TOML run file per argument, with the lattice and model in `[parameters]` and libraries in `input.lattice_library`/`input.model_library`; `lattice2xml LATTICE KEY=VALUE ...` is unchanged |
| `convert2xml`, `pyalps.convert2xml` | Load HDF5 results with `pyalps.loadSpectra`, `pyalps.loadEigenstateMeasurements` or `pyalps.loadMeasurements` |
| `alps::mcmpiadapter`, `alps::stop_callback(communicator, timelimit)` | `alps::mc::main` with `execution.chains` and `execution.time_limit` |
| `tutorials/09-code/07-mcmain-mcbase/mpi_pscan.cpp` and the `tutorials/10-ngs/3_mpi` and `4_mpi_pscan` programs | `tutorials/09-code/07-mcmain-mcbase/mpi.cpp`, one run file per scan point |
| `sparsediag`/`fulldiag`/`dmrg` XML job files and `--write-xml` output; `LATTICE_LIBRARY`/`MODEL_LIBRARY` parameters | TOML runs with `pyalps.run_io.execute(application, runs)`; `input.lattice_library`/`input.model_library` |
| XML-only averages: `sparsediag` "Ground State Energy"/"Energy Gap", `fulldiag` thermodynamics for a run's `T` or `beta` | Lowest energies from `pyalps.loadSpectra`; `fulldiag_evaluate --T value` (a `T` parameter of the run is its default) |
| `fulldiag_evaluate` on scheduler results (`.out.xml`, ALPS ≤ 3.0 `.out.h5`) | Rerun `fulldiag` with a TOML run; `pyalps.loadSpectra` still reads the old spectra |
| Several tasks of one XML job file in one program invocation | One TOML run per invocation; run a job manifest with `pyalps.run_io.execute(application, manifest, concurrency=N)` |
| `dmrg` `TEMP_DIRECTORY` parameter | `execution.temporary_directory` |
| `dmrg` `START_SWEEP`/`START_DIR`/`START_ITER` resumption from an earlier run's temporary files | Removed without replacement; rerun the calculation |
| Parapack `--threads-per-clone`, `--total-threads` | No thread-level replacement; run independent chains as MPI processes with `execution.chains` and `mpiexec -n` |
| Fortran F77 callbacks (`alps_init_`, ...), `ALPS_LONG`, scheduler input and XDR dumps | `bind(C)` callbacks taking a `type(c_ptr), value :: caller` handle, NUL-terminated names, `ALPS_INT64` and TOML runs; see the [Fortran examples](tutorials/00-examples/fortran/README.md#writing-an-application) |
| `src/apps/qmc/sse` (`dirloop_sse_v1`), `src/apps/qmc/sse2` sources | `dirloop_sse`; `sse2`'s optimized-ensemble workflow has no replacement |
| Top-level `import mpi` compatibility module | `import pyalps.mpi`; install the `mpi` extra for mpi4py |
| `applications/`, `tool/`, `test/`, `example/` | `src/apps/`, `src/tools/`, `tests/`, `tutorials/00-examples/` |
| Loose subsystem trees and transitional `src/alps/{common,runtime}/` | Semantic `src/alps/<module>/{include,src,tests}` ownership; see the [module map](src/alps/README.md#source-ownership) |
| `src/ietl/` | `src/alps/ietl/include/ietl/`; public `<ietl/...>` includes are unchanged |
| `<alps/scheduler/measurement_operators.h>` | `<alps/model/measurement_operators.h>` |
| `src/alps/config.h.in`, `src/alps/version.h.in`, transitional `src/alps/common/config/` | `cmake/config/` configuration templates |
| Build-tree `src/alps/` generated headers | `<build-dir>/generated/include/alps/`; consumers use exported CMake targets |
| XML definitions and stylesheets in `lib/xml/` | `src/alps/resources/` in the source tree; installed under `share/alps/xml/` |
| `tutorials/alpsize-*`, `tutorials/code-*`, `tutorials/ngs/` | `tutorials/08-alpsize/*`, `tutorials/09-code/*`, `tutorials/10-ngs/`; the Python export example is now in `python/pyalps/examples/ising/` |
| Root `ALPS_VERSION.txt` | `cmake/ALPS_VERSION.txt`, still shared by the SDK and Python package |

### Fixed

- Correct the Ising tutorials' Metropolis sign to match their ferromagnetic energy, accepting flips with probability `min(1, exp(-beta * delta_energy))`.
- Print real parameter values with the fewest digits that read back exactly, so TOML run-file couplings reach the lattice and model engines unrounded instead of at six significant digits.
- Keep model-library couplings and matrix elements exact. The expression engine reparses the text of partly evaluated terms and printed their numbers at six significant digits, so couplings and factors such as `sqrt(2)` entered Hamiltonians and measurement operators with relative errors near `1e-6`.
- Regroup every level in `pyalps.groupSets` for nested input; grouping already grouped data returned it unchanged, which broke the next-nearest-neighbor chain tutorial's critical-point analysis.
- Reclaim HDF5 variable-length string buffers and vector-attribute parent handles on conversion failure, preserving the conversion exception and partial-selection memory extents.
- Consolidate HDF5 writes with validation before replacement: reject malformed ranks, out-of-bounds or overflowing extents, and null buffers for nonempty transfers. Correct scalar/array attribute replacement and root attributes; restore serializer contexts when user hooks throw.
- Reject programmatic MaxEnt and CT-QMC input/output aliases after resolving paths against the preparation directory.
- Check all schema-declared input/output paths across TOML jobs before execution, including DMFT's final Green-function files and active text-output directories.
- Preserve empty run provenance in checkpoints, and reject NUL-containing parameter names and string values before overwriting stored parameters.
- Honor CT-HYB's resolved `VERBOSE` setting during progress reporting instead of reading an uninitialized flag.
- Prepare pinned toml++ headers for Ubuntu 22.04 CI, whose package archive does not supply the required parser.
- Report MaxEnt CLI help and input errors with normal exit codes instead of continuing into an empty input or aborting on an exception. Valid scientific runs are unchanged.
- Reject unsupported native C++ parameter checkpoint datatypes with the dataset path in the diagnostic, instead of silently substituting zero. Existing supported types and custom readers retain their decoding behavior; a failed parameter reload preserves the previous values.
- Generate XML plotting scripts compatible with Python 3 through the `alps-xml` CLI.
- Run the remaining tutorial notebooks under Python 3: the DMFT-02 analysis, the Japanese Python primer, and DMFT-08's headings stored as code cells. Notebooks and scripts no longer rely on invalid string escapes in plot labels, and the notebooks no longer refer to the retired `alpspython` launcher.
- Correct the Heisenberg tutorial's vector dot product and allow its vector implementation to compile without x86 SIMD support.
- Remove undefined behavior in XDR callbacks and integer-range boundary checks exposed by sanitizers.
