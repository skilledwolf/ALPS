# Module layout and ALPSCore reconciliation

ALPS sources are grouped by responsibility to prepare MaxEnt, HDF5 and typed params for ALPSCore reconciliation. Utilities, HDF5, params, Osiris, XML and command-line parsing have separate runtime libraries. Other source modules contribute to the existing `ALPS::alps` library or its shared compile interface. The first consolidation adapts the ALPSCore dictionary/value model in `params`; MaxEnt now reads TOML configurations through `run_config`. Scientific algorithms are preserved.

The [separate-process probes](../../tests/reconciliation/README.md) record the pinned ALPSCore reference and reproduce measured archive/params compatibility checks.

The HDF5 backend uses HDF5 2.x and HighFive 3.3.0. Complex values use compound
`r`/`i` members, Booleans use the `FALSE`/`TRUE` enum, scalars have rank zero, and
empty arrays retain explicit zero-length dimensions. The old marker datatypes,
trailing complex dimension and numeric/string cast engine are removed. Legacy
conversion belongs to the [standalone converter](../tools/hdf5/README.md), which
also upgrades explicitly versioned parameter checkpoints from v1 to v2. Generic
NULL dataspaces have lost their array shape; converting those or an application's
checkpoint layout still requires its scientific schema.

Use these buffer ownership, shape/type and error contracts with the new ALEA
maintainer. Start one measurement-application pilot alongside params
reconciliation. MaxEnt remains an
acceptance test for the foundations; it does not exercise ALEA migration.
`src/alps/alea/` still contains the legacy `Observable`/`ObservableSet`
implementation; this pass imports neither new ALEA nor a compatibility shim.
The pinned ALPSCore probes disable ALEA, so they do not validate ALEA compatibility.

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
belongs to utilities and solver execution. The shared archive registry remains
necessary for simultaneous handles to see replacement-file writes and for
read-only handles to share a later writer's context.

The measurement pilot in [checkpoint contracts](../../tests/pyalps/test_checkpoint_contracts.py)
uses deterministic scalar/vector sample streams. A checkpoint must preserve an
accumulator's unfinished bin so resumed statistics agree with an uninterrupted
stream. Aligned independent runs must preserve sample counts and means; their
merged errors are checked against an independent bin-level calculation. Results
must retain count, mean, error, variance, autocorrelation and bin data through
HDF5 reload, including scalar/vector shape and escaped observable names. These
tests constrain a future implementation without choosing between Core's
`accumulators` and separate newer `alea` APIs.

Uneven or incomplete legacy bin merging is an unresolved scientific gate. Direct
`MCScalarData` merging of samples 0–999 and 1000–2999 changes count from 3000 to
2976 on reanalysis and gives mean 1492.8333 instead of 1499.5. `mcobservable.merge`
of samples 0–516 and 517–1030 keeps count 1031 but gives mean 514 instead of 515.
These paths are unchanged from upstream master at `c22bfd701`; their outcomes
are not acceptance oracles for a replacement.
Require exact retained sample counts and means before admitting those merges;
errors for independent runs need not equal errors for a concatenated correlated
stream. Measurement collection loading also needs an explicit replacement-versus-
overlay contract. Existing result-format coverage does not settle these questions.

## Source ownership

Each module uses `include/`, `src/` and `tests/` where applicable. Public include spellings describe the API, independently of the physical owner: for example, NGS measurement headers live in `alea/include/alps/ngs/`, while typed parameters live in `params/include/alps/`.

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
| `alea/`, `accumulators/` | Observable/result facilities and accumulator implementations | `ALPS::alps` |
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

With application builds enabled, `alps-module-architecture.json` inventories 37 source owners, 509 public include spellings and 640 production files. The owners include `numeric_io`, `cli`, `plotting`, separate MaxEnt solver/executable owners, and nine [tool groups](../tools/README.md). Tool ownership includes historical inactive C++ sources without adding executable targets. These are ownership counts, not counts of independent libraries or passing tests. The earlier code checkpoint `f6f4501c0` had 24 owners, before the CLI and plotting modules were separated.

The foundation include cycle involving containers, HDF5, numerics, utilities and XML is removed. The current observed include graph retains the separate two-module cycle between `expression` and `legacy_parameters`; it still needs deliberate reconciliation. Dependency declarations constrain new include edges. Regenerate the report after changing module ownership or dependencies.

The report inventories 82 exact unresolved file/include pairs already present at baseline `0f7b995d5`: 80 in dormant accumulator code and two in the optional `USE_LATTICE_CONSTANT_2D` graph backend. Each exemption names its file, include and reason; they do not establish support for those inactive paths. Resolve or remove these dependencies deliberately rather than adding broad exclusions.

With `ALPS_BUILD_TESTING=ON`, CTest runs `module_architecture` and writes `<build-dir>/alps-module-architecture.json`; this requires a Python interpreter ≥ 3.11. Builds with testing disabled do not need Python for module configuration or manifest generation.

## Validation

`ALPS_BUILD_TESTING` controls module-local and central native tests. `ALPS_BUILD_APPLICATIONS` additionally controls MaxEnt and its regression. After building, focused existing labels include:

```sh
ctest --test-dir <build-dir> --output-on-failure -L '^(utility|hdf5|params|osiris|parser|cli|maxent)$'
```

Run the full native suite for shared-interface changes. The SDK consumers in `tests/cmake/consumer/` exercise aggregate and component links; integration tests for older inputs and observables keep `ALPS::alps`. Check shared/static consumers, installed-header ownership, SDK relocation and Python extension interoperability. The MaxEnt CLI contract checks TOML schema validation, help, and rejection before output creation. The MaxEnt regression also exercises early stop, callback exceptions and completion and cancellation semantics. These are validation requirements; this ownership map does not assert that every platform/configuration has passed them.
