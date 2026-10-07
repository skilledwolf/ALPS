# Testing ALPS

ALPS uses GoogleTest for C++ runtime tests, CTest for native execution, and
pytest for Python. Component tests live beside their implementation in
`src/alps/<component>/tests`; cross-component and installed-package checks live
here. The two language suites have independent entry points: passing CTest does
not mean that the Python package was tested.

## Local workflows

With the SDK prerequisites installed:

```sh
cmake --preset dev
cmake --build --preset dev --parallel 2
ctest --preset dev
ctest --preset dev -L hdf5
ctest --preset dev -R Matrix
```

`dev` omits applications. `default` includes application tests. `mpi` enables
MPI, `extensive` enables graph and HDF5 type-matrix tests, and `sanitizers`
instruments a Debug build with AddressSanitizer and UndefinedBehaviorSanitizer.
Build each preset before testing it. To exercise actual ranks:

```sh
cmake --preset mpi
cmake --build --preset mpi --parallel 2
ctest --preset mpi -L '^mpi$'
```

All test presets reject an empty selection. An unbuilt GoogleTest executable is
reported as `NOT_BUILT`; a built executable with zero discovered cases fails
discovery, even when other executables still contain tests.
GoogleTest is a test-only dependency: CMake uses an installed GTest >= 1.14 or
fetches the checksum-pinned 1.18.0 release. An offline source checkout can be
selected with `-DFETCHCONTENT_SOURCE_DIR_GOOGLETEST=/path/to/googletest`.
`ALPS_BUILD_TESTING=OFF` neither finds nor fetches GoogleTest, and it is never
installed as part of the SDK.
Tutorials using only the CMake process runner do not acquire GoogleTest either.

Standalone numerical tutorials need no SDK or wheel:

```sh
python -m pip install numpy pytest
python -m pytest tests/tutorials -q
```

After installing the SDK and pyalps as described in `CONTRIBUTING.md`:

```sh
python -m pytest tests/pyalps -q -rs
python -m pytest tests/cmake -q -rs
python -m pytest tests/ci tests/packaging -q
```

`tests/cmake` requires `ALPS_DIR`; it validates installed consumers, relocatability,
build configuration and the testing harness. `ALPS_TEST_CMAKE_ARGS` supplies a
JSON array of dependency/toolchain arguments to its temporary projects. Keep
the installed-SDK and Python checks pointed at the candidate build. Distribution
CI additionally exercises repaired wheels and independently unpacked sdists.

## Writing a C++ test

Create an ordinary executable, link only the component being tested, then call
`alps_add_gtest(target COMPONENT hdf5 LABELS compatibility)` from its test
CMakeLists. Each GoogleTest case is discoverable through CTest; names include
the component and executable to avoid collisions across generated type tests.
Use `TIMEOUT` for a justified nondefault limit and `THREADS` for an OpenMP test.
Labels describe both ownership and purpose (`compatibility`, `integration`,
`scientific`, `extensive`, `mpi`). GoogleTest execution alone does not imply a
test is a scientific validation.

Use `alps::testing::TemporaryDirectory` for files. Cases in one executable may
run concurrently in different processes, so an executable-level working
directory is not enough to isolate filenames. Tests must not depend on ordering
or prior test output. GoogleTest assertion failures remain active in Release
builds; ordinary C `assert()` is unsuitable for test expectations.

MPI tests use `alps_add_mpi_gtest` and run as complete rank topologies. The shared
main initializes MPI, combines rank exit statuses, and writes distinct XML
reports. Every rank must enter matching collectives even after nonfatal assertion
failures. Use CTest timeouts to bound deadlocks. The rank count is also declared
as the CTest processor requirement; OpenMP and BLAS threads are bounded separately.

## Choosing assertions and reference data

Prefer direct value assertions to printed numbers. For floating-point results,
state the absolute/relative tolerance and its numerical or physical basis.
Historical decimal fixtures migrated to assertions retain their original
rounding bounds; those comparisons establish compatibility, not independent
scientific correctness. Avoid broadening a tolerance just to pass a new platform.
Floating-point equality and tolerance assertions reject NaNs; test intended
nonfinite behavior explicitly.

Combine reference values with analytic solutions, independent calculations and
invariants such as normalization or conservation. A writer/reader round trip
can hide a shared bug: retain historical files and raw-schema assertions too.
Fixture changes require provenance and a behavioral explanation in the PR.

Seed stochastic regression tests and include the seed in diagnostics. These
reproduce a particular trajectory; they do not validate uncertainty estimation
in general. New statistical validation should document the ensemble,
autocorrelation treatment, thresholds and combined false-failure probability.
Do not retry failed numerical tests until they happen to pass.

Textual serialization contracts can use `StreamFixture`: it supplies an input
fixture and compares captured output explicitly inside GoogleTest. Existing
model/lattice transcripts are retained as compatibility evidence. They should
be supplemented with semantic assertions when changing the underlying behavior.
The remaining CMake process runner is for executable/CLI/tutorial contracts.

## CI policy and evidence

Every PR runs the standalone NumPy suite at its floor and current version, plus
CI/packaging helper checks. Native and packaging jobs use conservative path
classification: standalone optical-lattice tutorial/prose-only changes can take
the fast path; ordinary implementation changes use representative configurations;
shared headers, dependencies, CMake/CI changes and unknown paths select the full
matrix. Missing diff metadata also selects the full matrix.

Master pushes, release tags, weekly schedules and manual runs select the full
supported matrix. Native jobs include real MPI execution when enabled. A
dedicated sanitizer job exercises the native suite. Stable aggregate checks
(`Source CI` and `Packaging CI`) require all selected jobs to succeed; these are
the checks to select in repository branch protection. JUnit reports, CTest logs,
rank reports and
configuration diagnostics are retained for inspection.

Benchmarks remain separate from correctness gates. Existing timing programs
are not claims of performance coverage. Introduce baselines only with a defined
workload and controlled measurement environment; ordinary shared-runner wall
times are diagnostic, not reliable regression thresholds.

## Migration accounting

The migration preserves existing executable identities where possible and
expands Boost.Test typed cases into individually discovered GoogleTest cases.
Counts before and after migration are therefore not directly comparable.
Component migration notes record retained scenarios, reference provenance and
historical behavior needing a consolidation decision:

| Area | Contract and inventory |
| --- | --- |
| Numerics, containers, alea, legacy parameters, parapack | [Scenario and reference accounting](migration.md), including numerical quirks preserved for an explicit consolidation decision. |
| HDF5 | [Schema, fixtures, generated types and resource opt-ins](../src/alps/hdf5/tests/README.md). |
| Params | [Typed access, persistence and adapter contracts](../src/alps/params/tests/README.md). |
| Graph | [Canonicalization, embedding and extensive-test inventory](../src/alps/graph/tests/README.md). |
| MPI | [Rank topologies, assertions and remaining manual drivers](mpi.md). |
| Model, lattice, XML, Osiris | Existing transcript oracles remain exact GoogleTest compatibility checks; Release-disabled assertions are replaced and dump files are isolated. |
| Utilities, CLI, portability and random | Direct API, wire-format, seeded sampling and platform assertions replace ad hoc mains or printed smoke output. |
| DMFT, DMRG, MaxEnt | Native regression tests use GoogleTest; executable behavior stays under CTest. MaxEnt's Python scientific oracle is shared with its binding tests. |
| Monte Carlo runners | The formerly unregistered scalar/vector examples check completion, stopping, collection and HDF5 persistence. Test counters are initialized; MPI scheduling is deterministic. |

The obsolete experimental accumulator prototypes were removed in #163.
Benchmarks and remaining unregistered sources do not count toward executed coverage.

Use named skips for unavailable optional capabilities or unsupported type
combinations; skipped cases do not establish coverage. Do not add permanently
disabled tests. Any necessary quarantine must name the
reason, tracked issue and owner; CI must continue reporting the missing coverage.
