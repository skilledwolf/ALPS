# Module layout and ALPSCore reconciliation

ALPS sources are grouped by responsibility to prepare MaxEnt, HDF5 and typed params for ALPSCore reconciliation. Utilities, HDF5, params, Osiris, XML and command-line parsing have separate runtime libraries. Other source modules contribute to the existing `ALPS::alps` library or its shared compile interface. The first consolidation adapts the ALPSCore dictionary/value model in `params`; MaxEnt now reads TOML configurations through `run_config`. Scientific algorithms are preserved.

The [separate-process probes](../../tests/reconciliation/README.md) record the pinned ALPSCore reference and reproduce measured archive/params compatibility checks.

The HDF5 backend uses HDF5 2.x and HighFive 3.3.0. Complex values use compound
`r`/`i` members, Booleans use the `FALSE`/`TRUE` enum, scalars have rank zero, and
empty arrays retain explicit zero-length dimensions. The old marker datatypes,
trailing complex dimension and numeric/string cast engine are removed. Legacy
conversion belongs to the [standalone converter](../tools/hdf5/README.md), which
also migrates selected released parameter and ALEA schemas. Generic
NULL dataspaces have lost their array shape; converting those or an application's
checkpoint layout still requires its scientific schema.

The Eigen-based modern ALEA core is independently linkable as `ALPS::statistics`.
Its [serialization contracts](alea/modern-alea.md) and the accumulator-only Ising
pilot cover results, complex covariance and actual batch-accumulator continuation.
The thin HDF5 adapter shares the canonical primitive mappings. Legacy
`Observable`/`ObservableSet` APIs remain in `ALPS::alps` while active clients migrate.
Standalone CT-INT, CT-HYB and Hirsch-Fye use modern ALEA, pooling raw independent
replicas before signed analysis and publishing canonical results atomically.
CT-INT includes general multiband density interactions and their full density
moments; DMFT uses only external solver processes. Other legacy measurement
producers still need migration before their statistical interfaces can be deleted.
The NGS `mcbase` framework and its MPI adapter use the same native batches,
checkpoints and result codec. Python measurement maps share native accumulator
handles; collected result dictionaries own snapshots. The old NGS measurement
facade, feature-stack API and unused scheduler prototypes are removed.
Historical reconciliation probes disabled ALEA; they do not validate this import.

## Consolidation acceptance contracts

The current params layer adopts Core's owning dictionary/value model with deliberate
ALPS value semantics. Signed and unsigned integers use 64-bit storage; real and
complex values use double precision, alongside Boolean, string and homogeneous
vector alternatives. Conversions check range and precision. Boolean-to-integer,
string-to-number, real-to-integer and scalar-to-vector access are rejected;
exact integer-to-real and elementwise vector conversions are supported. Copies
own their values. Nonconst missing lookup inserts an unset entry, while const
lookup throws; `exists()` excludes unset entries. Iteration is sorted and erasing
an absent name is harmless. The [native contracts](params/tests/contract.cpp)
and [Python conversion contracts](../../tests/pyalps/test_conversion_contracts.py)
are acceptance tests for these choices.

Schemas, defaults and provenance belong to `run_config`, rather than the
dictionary. `alps.params.v2` is the canonical parameter checkpoint for this
branch; it deliberately rejects older ALPS/Core parameter encodings. Generic
archive payload interchange is a separate contract. The reconciliation runner
checks both providers' own extended checkpoints when both are supplied, and can
check the installed ALPS SDK alone. A same-provider report does not demonstrate
Core interoperability, and the historical comparison is not the current baseline.

Archive reads copy values into caller-owned storage and must release HDF5-allocated
variable-length buffers on success and assignment failure. Partial selections
must use the matching memory extent for cleanup. Cleanup must preserve the
original exception. The existing `hdf5_valgrind` regression exercises
scalar/vector datasets and attributes, including partial selections.
HighFive owns files, groups, datasets, datatypes, dataspaces and attributes; the
public overloads enter one locked implementation. The
`hdf5_read` contract also covers all native numeric source types, fixed-width
strings, selections and rejected reads. HDF5 performs numeric conversions;
cross-family string/numeric and Boolean/numeric conversions are rejected. Direct
HighFive and h5py fixtures test interoperability independently of ALPS writers.

Writes share object creation/replacement, layout selection and transfer code.
Rank, bounds, extent products and nonempty buffer pointers are checked before
replacing stored data. Attribute replacement respects datatype and shape, including
scalar-to-array changes; partial attribute transfers are rejected. Internal HDF5
handles have RAII ownership. Object serializers restore their previous context
on both return and exception; archive copies share file ownership, while copy
assignment is disabled because overwriting the context pointer bypasses that ownership.

Opening an archive leaves process signal handlers alone. Termination polling
belongs to utilities and solver execution. Independent opens have independent
permissions; copied views share ownership and close together.

The [checkpoint contracts](../../tests/pyalps/test_checkpoint_contracts.py)
cover complete scalar/vector batch continuation, RNG state, escaped measurement
names and failed-load preservation. NGS base loads stage parameters, measurements
and RNG together, require already registered measurements and validate their
component dimensions. Application subclasses validate their own spin/progress
state before calling the base load. Virtual archive hooks take an archive reference;
the former by-value tutorial hooks did not override the base and skipped that state.

MPI collection checks result names on every rank and includes empty ranks in
raw reduction. Independent partial bins retain their original weights; pooling
never pairs unrelated replica bins. Scheduled failure consensus propagates local
sampling/callback errors before progress collectives. The native Ising and
Heisenberg contracts verify exact restart and independently calculated Boltzmann
energies; they do not use legacy statistics as a correctness oracle.

Direct legacy `MCScalarData` still serves other applications and has unresolved
uneven-bin merging: samples 0–999 and 1000–2999 change count from 3000 to 2976
on reanalysis and give mean 1492.8333 instead of 1499.5. This upstream behavior
is a migration gate for those remaining callers. Native result reduction retains
exact sample counts, sums and weights. Released recoverable linear-bin histories
can be converted offline to native batch analysis; missing joint covariance or
checkpoint cursors cannot be reconstructed from summaries.

## Source ownership

Each module uses `include/`, `src/` and `tests/` where applicable. Public include spellings describe the API, independently of the physical owner: for example, the simulation framework lives in `mc/include/alps/`, while typed parameters live in `params/include/alps/`.

| Source module under `src/alps/` | Responsibility | Binary or compile owner |
| --- | --- | --- |
| `utilities/` | Utility functions, general helpers, type traits, termination polling and NGS configuration helpers | `ALPS::utilities` |
| `containers/` | Fixed-capacity containers, ALPS multi-array storage and Boost serialization | `ALPS::containers` |
| `numerics/` | Numerical helpers, array mathematics and matrix/vector interfaces | `ALPS::numerics` |
| `numeric_io/` | HDF5 adapters for numerical matrices and vectors | `ALPS::numeric_io` |
| `ietl/` | Iterative eigensolver headers under `include/ietl/` | `ALPS::headers` |
| `hdf5/` | Archive API, container adapters and explicit checkpoint publication | `ALPS::hdf5` |
| `params/` | ALPSCore-derived owning values, checked lookup/conversion and versioned HDF5 checkpoints | `ALPS::params`; `adapters/` contributes to `ALPS::alps` |
| `run_config/` | TOML run files, application schemas, defaults, validation and provenance | `ALPS::run_config` |
| `osiris/` | Dump serialization, process/communication state and XDR implementation | `ALPS::osiris` |
| `xml/` | XML parsing, handlers, attributes and output streams | `ALPS::xml` |
| `cli/` | Existing `mcoptions` and `parseargs` command-line grammars | `ALPS::cli` |
| `plotting/` | `<alps/plot.h>` output helpers combining XML and older parameters | `ALPS::headers` |
| `legacy_parameters/`, `expression/` | Older `alps::Parameters` and expression evaluation | `ALPS::alps` |
| `graph/`, `lattice/`, `model/` | Graph helpers, lattice definitions and physical models | `ALPS::headers`, `ALPS::alps` |
| `random/` | Random generators and their factories | `ALPS::alps` |
| `alea/` | Legacy observables and modern Eigen-based statistical estimators | `ALPS::alps`, `ALPS::statistics` |
| `mc/`, `scheduler/`, `parapack/` | Simulation API, execution and scheduling | `ALPS::alps` |
| `fortran/` | C++ bridge with public headers in `include/alps/fortran/` | `ALPS::fortran` |
| `solvers/` | Shared `<alps/solvers.hpp>` declarations for MaxEnt and CT-QMC | `ALPS::solver_headers` |
| `resources/` | XML definitions and stylesheets | Installed data component |

The broad `common/` and `runtime/` source groups are removed. ALPS and IETL configuration templates live in `cmake/config/`; generated public headers live under `<build-dir>/generated/include/alps/`. `ALPS::configuration` carries the shared configuration headers and compile requirements. MaxEnt remains in `src/apps/maxent/{src,cli,tests}`. Subsystem tests follow their source owner; cross-module integration, SDK, Python, CLI, packaging and reconciliation checks stay under root `tests/`. Relocation preserves registered test names and leaves inactive fixtures inactive.

All first-party public headers have explicit CMake `HEADERS` file sets. These declare build include roots and preserve installed `<alps/...>` and `<ietl/...>` paths; consumers use exported targets rather than broad source-tree include roots. Public template definitions remain installed, while private source files, tests and physical `include/` nesting do not enter the SDK.

## Library boundaries

`ALPS::configuration`, `ALPS::containers`, `ALPS::numerics`, `ALPS::numeric_io` and `ALPS::solver_headers` are exported interface targets, not new binary libraries. Foundation targets own their public header sets and declare their actual dependencies; they do not inherit `ALPS::headers`. The aggregate compile interface instead collects these foundations and the remaining simulation headers.

`ALPS::containers` provides storage without the mathematical umbrella. Array functions/operators, `<alps/multi_array.hpp>`, `<alps/functional.h>` and `<alps/utility/numeric_cast.hpp>` now belong to numerics, with their implementations and public include names preserved. Utility array resizing, abbreviated printing and MPI helpers include only the container facilities they use. `ALPS::numerics` adds utilities, containers and BLAS/LAPACK without HDF5 or XML.

`ALPS::numeric_io` combines numerics and HDF5. Include `<alps/hdf5/matrix.hpp>` or `<alps/hdf5/numeric_vector.hpp>` explicitly and link `ALPS::numeric_io` when using these archive adapters; `<alps/numeric/matrix.hpp>` no longer includes the matrix archive adapter. Diagonal/deprecated matrix/vector HDF5 `save`/`load` members and matrix `write_xml`/XML insertion overloads, unused in this repository, are removed without replacement adapters. This pass therefore includes API removals as well as ownership moves; it changes no numerical algorithms and imports no ALPSCore implementation.

`ALPS::utilities` owns utility symbols and the NGS termination-signal queue. It links Boost.Filesystem and platform threads without the simulation runtime, HDF5 or BLAS/LAPACK.

`ALPS::hdf5` owns archive symbols, exception exports and file ownership. It uses HighFive headers privately and links utilities, HDF5 and platform threads. Boost.Filesystem is supplied by utilities. Installed consumers do not need HighFive headers or its CMake package. Independent opens keep their own permissions; copies are views closed together by explicit `close()`. Modes are exactly `r`, `a` and `w`; `save_checkpoint` publishes a complete file only after serialization and checked close. Legacy formats belong to the standalone converter. Python ordinary IO uses h5py, temporarily transferring file ownership for native scientific serialization. A parallel HDF5 provider can bring its own MPI dependency.

`ALPS::params` owns ALPSCore-derived dictionary/value storage and explicit `alps.params.v2` checkpoints. Native datatype and scalar/vector rank identify values; indexed name/value entries preserve parameter names without duplicate type tags. It links HDF5 and Boost.Serialization; MPI builds also use MPI and Boost.MPI. Python values are eagerly copied into native storage. The file constructor, XML reader, proxies and Python `paramvalue_source` interface are removed. The remaining `params/adapters/` function converts typed scalars to the older model/lattice `Parameters` API for live internal callers; it is compiled into `ALPS::alps`.

`ALPS::run_config` applies separate application TOML schemas, defaults, type/range checks, relative-path resolution and run provenance. It uses toml++ headers privately in one translation unit, including when the package manager provides a compiled toml++ library. Parsing and application orchestration do not belong to the dictionary.

`ALPS::osiris` owns dump/process APIs, XDR symbols and the communication state used by `comm_init()` and `is_master()`. It links Boost.Serialization/Filesystem and, when enabled, MPI. This gives communication state one owner and permits MaxEnt to preserve existing diagnostic gating without linking the simulation runtime.

`ALPS::xml` owns XML parsing, attributes, handlers, output streams and stylesheet lookup, with Boost.Filesystem and Boost.Regex dependencies. The separate `plotting/` header owner keeps `<alps/plot.h>` and its older `Parameters` dependency outside the XML component.

`ALPS::cli` owns the existing `mcoptions` and `parseargs` implementations, linking utilities and Boost.ProgramOptions. Their installed headers remain `<alps/ngs/mcoptions.hpp>` and `<alps/parseargs.hpp>`. The two existing option grammars, defaults, filename rules and error behavior are preserved; these parsers do not read parameter files.

`ALPS::maxent` links run configuration, params, HDF5, utilities, Osiris and its Boost/numerical providers. Numerical algorithms are preserved. Its entry point takes resolved parameters, a separate scientific-data object and explicit output/execution settings; its return value reports whether a result was produced. The CLI reads a TOML run file and supports validation before output creation. Neither the solver nor executable links `ALPS::alps`. Public declarations are in `<alps/maxent.hpp>` and `<alps/solvers.hpp>`.

`ALPS::alps` links the extracted runtime components publicly. These libraries follow `BUILD_SHARED_LIBS`, with component-specific generated export headers. Python extensions require shared runtime libraries and package one copy of every component in `ALPS_RUNTIME_TARGETS`. Rebuild downstream binaries after the XML and CLI extractions, as after the earlier library splits.

Physical ownership does not imply independent linkability for every module. `ALPS::headers` still exposes the aggregate compile interface for the simulation modules, including the remaining expression/older-parameters cycle. The foundation interfaces have narrower dependencies, but package discovery still checks the full SDK dependency set. Shared build policy remains in the root CMake files and `cmake/`.

## Architecture checks

CMake generates `alps-module-manifest.json` from module declarations, public-header file sets and source lists. The [architecture checker](../../.github/scripts/check_module_architecture.py) checks production-file ownership, public include ownership, declared include dependencies and public/private boundaries. It reports observed cycles, nonliteral includes and documented unresolved includes for review. It does not replace compilation, link-dependency checks or scientific validation.

After configuring the build, run:

```sh
python .github/scripts/check_module_architecture.py \
  --manifest _build/default/alps-module-manifest.json \
  --write-report _build/default/alps-module-report.json
```

Update the owning module's CMake declarations when adding files or dependencies; the manifest is generated rather than hand-maintained. A declared existing cycle is visible architectural debt, not evidence of independent components.

### Recorded architectural debt

The architecture manifest records current source ownership, public includes and
production files. Tool ownership includes historical inactive C++ sources;
ownership does not establish independent linkability or passing tests.

The foundation include cycle involving containers, HDF5, numerics, utilities and XML is removed. The current observed include graph retains the separate two-module cycle between `expression` and `legacy_parameters`; it still needs deliberate reconciliation. Dependency declarations constrain new include edges. Regenerate the report after changing module ownership or dependencies.

Only the two exact unresolved includes in the optional `USE_LATTICE_CONSTANT_2D`
graph backend remain exempted. The dormant accumulator implementation and its
80 obsolete include exemptions are removed.

With `ALPS_BUILD_TESTING=ON`, CTest runs `module_architecture` and writes `<build-dir>/alps-module-architecture.json`; this requires a Python interpreter ≥ 3.11. Builds with testing disabled do not need Python for module configuration or manifest generation.

## Validation

`ALPS_BUILD_TESTING` controls module-local and central native tests. `ALPS_BUILD_APPLICATIONS` additionally controls MaxEnt and its regression. After building, focused existing labels include:

```sh
ctest --test-dir <build-dir> --output-on-failure -L '^(utility|hdf5|params|osiris|parser|cli|maxent)$'
```

Run the full native suite for shared-interface changes. The SDK consumers in `tests/cmake/consumer/` exercise aggregate and component links; integration tests for older inputs and observables keep `ALPS::alps`. Check shared/static consumers, installed-header ownership, SDK relocation and Python extension interoperability. The MaxEnt CLI contract checks TOML schema validation, help, and rejection before output creation. The MaxEnt regression also exercises early stop, callback exceptions and completion and cancellation semantics. These are validation requirements; this ownership map does not assert that every platform/configuration has passed them.
