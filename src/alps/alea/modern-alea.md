# Modern ALEA statistical core

`ALPS::statistics` contains the Eigen-based ALEA statistical core imported from
ALPSCore revision `7146b9e1f017938a94e5dae35d88467cc5ba7969`. The source and
small `common::ndview` / serialization support headers retain their original
ALPS Collaboration copyright notices and MIT licensing. Core's MPI and stream
codec plugins and package build system are not imported. Existing legacy
ALEA APIs remain with `ALPS::alps` while active clients migrate.

## ALPSCore consolidation decisions

Native ALEA is the shared statistical implementation. The older ALPSCore
`accumulators` package is outside this consolidation: its wrapper API and
second estimator stack will not be imported. ALPS uses its own HDF5 codec and
Boost.MPI reduction bridge; the upstream codec plugins and build system would
duplicate those facilities.

`alps::mc` owns TOML application execution, independent chains, collective
workers, atomic publication and restart. `mcbase` remains a simulation API
used by those applications and the teaching examples; its `run`/`mcmpiadapter`
entry points also remain available for custom embedding. These are deliberate
ALPS interfaces, not a promise of drop-in ALPSCore Monte Carlo API compatibility.
Tutorial front ends still need consolidation onto the shared TOML runner.

### Changes from the pinned ALEA source

Keep this ledger current when updating the import. Upstream provenance is the
revision above; this is a maintained integration, not an unmodified vendor copy.
The original [license](LICENSE.TXT), [copyright](COPYRIGHT.TXT) and
[acknowledgment](ACKNOWLEDGE.TXT) notices accompany source and SDK distributions;
the combined MIT notice also ships in Python packages.

| Area | Change and reason | Regression coverage |
| --- | --- | --- |
| Moment accumulation and merging | Centered Chan/Welford sums replace cancellation-prone raw subtraction; real, circular and elliptic complex strategies retain their distinct algebra. Squared weights use floating arithmetic to avoid upstream `uint64` overflow. | `tests/statistics.cpp`, `tests/mpi.cpp` |
| Independent runs and MPI | Reduction stages owned snapshots instead of mutating const inputs. Empty runs are neutral; all ranks validate schemas before collectives. Batch histories retain independent-run weights. `reduce` is an ALPS bridge API. | `tests/mpi.cpp`, `tests/statistics.cpp` |
| Nonlinear propagation | Unequal/empty batches are handled explicitly; signed ratio covariance retains cross terms; the numerical Jacobian uses central differences. Declared but unimplemented bootstrap/sampling tags are removed. | `tests/statistics.cpp`, Python ALEA tests |
| Mean tests | Correct effective weights and dimensional degrees of freedom, use Boost's F distribution, and expose explicit one-/two-sample `t2_test` overloads. | `tests/statistics.cpp` |
| Results and restart | Native result codecs validate before replacement; kinds 6–10 add full accumulator checkpoints, including partial batches and hierarchy state. Loading a result into an autocorrelation accumulator is rejected: upstream's unfinished overload could not restore that state. | `tests/checkpoint_contracts.C`, `tests/migration.cpp` |
| Public integration | Native Python bindings, family joins, real-component transforms and diagnostics expose the same statistical core. Serialization adapters keep HDF5 outside `ALPS::statistics`. | Python ALEA tests, installed SDK consumers |

Use `<alps/alea.hpp>` for estimators and `<alps/alea/hdf5.hpp>` for the thin
serializer bridge to the canonical native HDF5 mappings; users of that adapter
also link `ALPS::hdf5`. The statistical core itself has no HDF5 dependency. Eigen matrices use
physical `[columns, rows]` axes; elliptic complex covariance retains its real
2×2 operator axes. Reads never create groups, and failed result loads preserve
the destination. Noncontiguous Eigen serialization evaluates once.

Each result group has unsigned scalar `@version=1` and `@kind`: mean 1,
variance 2, covariance 3, autocorrelation 4, batch 5. Existing Core scientific
field names, counts, squared weights and full covariance are retained. Derived
error fields remain available for inspection and are validated when reading.

`<alps/alea/checkpoint.hpp>` supplies actual resumable `batch_acc<T>`
`serialize`/`deserialize` overloads, kind 6. Checkpoints retain sums, per-batch
counts, offsets, base batch size and the complete merge cursor; loads stage and
validate that state before replacement. Only the cursor's level and position
are stored; its other fields are derived. The codec replays native batch merges
from the sample count to validate every count, offset and cursor together;
empty bins must have zero sums. Samples are unit observations; result bins retain
their accumulated sample weights. Merging unrelated time series into a
resumable batch accumulator is unsupported; variance and covariance estimators
retain weighted result merging.

The same header supplies resumable codecs for mean (kind 7), variance (8),
covariance (9), and autocorrelation (10), including real, circular complex and
elliptic complex moments. Moment checkpoints store centered sums without a
round trip through normalized results and retain the unfinished batch. Counts
include partial samples. Autocorrelation checkpoints retain each level's partial
state and validate the hierarchy against the total sample count, base batch
size and granularity. Restart therefore reproduces uninterrupted statistics.
Loading stages the complete state before replacement. Inserting an independent
result into a live autocorrelation accumulator remains unsupported: a result
cannot reconstruct a missing partial hierarchy.

Python exposes native `Mean`, `Variance`, `Covariance`, `Autocorrelation` and
`Batch` accumulator/result pairs, with `Complex` variants and `EllipticVariance`
and `EllipticCovariance`. Samples are scalar or 1D numeric arrays; means and
errors are vectors, including one-component results. Elliptic uncertainty arrays
retain explicit trailing `[real/imag, real/imag]` axes. Arrays own their data.
`result()` takes a snapshot. `save(archive,path)`, `load(archive,path)` and
`Class.read(archive,path)` use the native result or full checkpoint codec.
Batch accumulators expose offsets; autocorrelation results expose the hierarchy,
level estimates, `tau` and whether a coarse level is available.

`alea.merge(results)` pools independent runs with their original weights;
`result.join(other)` concatenates components, preserving aligned batch evidence
or assuming independence between summary results. Covariance joins retain
squared weights. `result.transform(function, output_size, method, dx)` calls the
native propagation code: `none`, `linear`, or `jackknife` where the evidence
supports it. Batch transforms default to jackknife; other uncertain estimates
default to linear propagation. A callback receives and returns 1D arrays.
`test_mean(expected)` accepts a reference vector or another result of the same
type and returns the native test statistic and probabilities. Complex batches
and elliptic full covariance can expose `real_components()` for arbitrary
real/imaginary transforms, retaining the joint evidence. Circular complex
summary covariance cannot reconstruct this information and rejects conversion.
A single elliptic variance component also supports this conversion;
`ratio_real_imag` handles componentwise signed ratios directly.

`pyalps.loadMeasurements` reads all five native result kinds. Each returned
`DataSet.native_result` retains the full native estimator for subsequent analysis;
`y` remains the convenient plotting projection. Elliptic plotting errors use the
circular magnitude while native results retain the exact 2x2 uncertainty blocks.

The NGS `mcbase` framework stores shared native accumulator handles in a standard
variant, allowing each measurement to choose any estimator family or supported
complex strategy. Python assigns accumulator objects directly to
`sim.measurements[name]`; a retrieved handle stays valid after replacement or
erasure. Collected dictionaries contain the concrete native result types. C++
can access `measurement<Accumulator>(name)` (real batches by default), inspect the
standard result variant, or request a homogeneous collection with
`collect_results_as<Result>()`. Spin engines explicitly select batch evidence
for their joint jackknife analysis. Base checkpoint loading stages parameters,
measurements and RNG and validates every registered name, type and dimension
before replacement. The existing native kind, datatype and shape identify the
stored estimator; there is no second type-tag protocol. Subclasses override
archive-reference hooks for their application state.


Variance and covariance accumulate centered weighted moments using Chan/Welford
updates. Merging never modifies its input result. MPI reductions sum local
centered moments and combine the gathered run means and weights, avoiding raw
second-moment cancellation. This applies to real, circular complex and elliptic
complex estimators, including the variance hierarchy used for autocorrelation.
The extra reduction storage is linear in ranks times components, even when the
estimator holds a full covariance matrix. Constant unit streams have zero error;
empty and single-observation streams still have unavailable variance estimates.

Independent-run result reduction retains every batch and its weight, including
unfinished batches. Runs with different slot counts use disjoint blocks padded
with empty slots; their bins are never summed together. Autocorrelation reduction
retains only levels present in every nonempty run; empty runs are neutral. Every
retained level contains all samples. Failed reductions preserve the original
result. Reduction does not create a resumable combined time series.

The optional `<alps/alea/mpi.hpp>` reducer takes a borrowed `MPI_Comm`; clients
link `MPI::MPI_CXX` alongside `ALPS::statistics`. The statistics library itself
remains MPI-free. Construction and reductions are collective and must occur in
the same order on every rank. Only the chosen root retains the combined result.
Sample counts use unsigned 64-bit MPI arithmetic.
`mcmpiadapter` checks matching result requests and estimator types before reducing every entry,
including empty local accumulators. Only the root receives collected results.
Local sampling and callback failures reach every rank at scheduled checks.
Custom reducers must implement `reduce(view<uint64_t>)`; rebuild downstream
binaries after this interface change. The combined sample count must fit in
`uint64_t`.

Signed estimates use joint batches of the signed numerator and sign. Apply the
existing binary transformer to their ratio, with linear covariance propagation
or weighted jackknife propagation; separate scalar results lose the covariance
needed for the error. The transformer defines its domain and should reject a
zero denominator. Jackknife requires two occupied bins and skips empty slots.
Its result is an analysis result, not an accumulator checkpoint.

For large collections of componentwise signed measurements, an elliptic
`var_acc<std::complex<double>, elliptic_var>` can store each signed numerator
as the real part and its sign as the imaginary part. `ratio_real_imag(result)`
returns a real variance result using their retained 2×2 covariance. Storage
is linear in the number of components; these errors assume independent samples
and do not estimate autocorrelation or covariance between components. Empty
and single-observation results retain unavailable errors; zero average signs
and invalid ratio domains fail without changing the source result.
The heterogeneous `result` facade also exposes `count2()` and `stderror<T>()`;
mean-only results cannot provide either estimate.

`test_mean` evaluates Hotelling's mean comparison after diagonalizing the full
covariance. The two-sample test assumes independent Gaussian samples with equal
population covariance; it pools covariance with `n1+n2-2` degrees of freedom.
Effective observation counts remain floating point and provide an approximation
for weighted bins or correlated samples. See the
[one- and two-sample formulas](https://online.stat.psu.edu/stat505/Lesson07).
Complex batches and elliptic results are expanded into joint real/imaginary
vectors. Circular-only results lack the covariance needed for this test and
are rejected; componentwise variance results are supported only for one scalar
or one elliptic complex value. Deterministic equal/mismatched directions produce probabilities
one/zero; invalid covariance and insufficient observations are rejected.
The low-level `t2_test(diff, variance_of_mean, covariance_dof, atol)` takes
diagonalized real data explicitly, without encoding degrees of freedom in a
synthetic result. F-distribution tails use Boost.Math.

The standalone CT-INT solver stores sign-weighted physical measurements as
joint native batches. It reduces those batches across MPI ranks before taking
ratios, keeping independent partial bins and their counts. Unsigned diagnostics
use ordinary batches. Its published analysis results use the same kind-5 codec;
`pyalps.loadMeasurements` reads their means and errors through the native reader.
The CT-INT driver publishes the result, Green functions and run configuration
together through the existing checked HDF5 publication helper. These files are
analysis outputs; CT-INT does not expose a solver restart interface.
The same density-density kernel serves two-flavor and general multiband runs.
It retains the full flavor-pair density moments, with occupation idempotence on
the diagonal, for general interaction-dependent Fourier tails in both measurement
modes. Uniform unordered interacting-pair proposals preserve detailed balance
for sparse matrices; isolated flavors and zero interactions need no retry loop.

CT-HYB uses the same native collection and publication path. Time measurements
pair their averaged numerators with the sign averaged over `N_MEAS` updates;
final-configuration measurements use that configuration's sign. Sign-weighted
bounded measurements use joint batches and weighted jackknife; Sign and order
diagnostics use ordinary batches. G2/H2 tensors use the componentwise ratio
above. They publish kind-5 and kind-2 results, respectively,
and `pyalps.loadMeasurements` retains both families' means and errors. Physical
G/F endpoints are formed before accumulation, so derived means, errors and
covariance share the same convention. These outputs do not contain the segment
configuration needed to restart CT-HYB.

Hirsch-Fye also uses native joint batches and raw MPI collection. Physical
Green samples include the beta endpoint before signed analysis, and warm-up
excludes the scheduler's former transition/reset sweep. Its canonical kind-5
results and derived time/frequency Green functions share one atomic output
with typed parameters and provenance. The auxiliary-spin chain has no solver
restart interface.

`pyalps.hdf5.save_checkpoint(filename, callback)` calls the existing native
publication helper. The callback receives a `NativeArchive`; retained callback
views close before publication, and a failed save preserves the previous file.
The pure Python Ising tutorial uses this path and checks complete restart against
an uninterrupted spin stream. Both Python and C++ pilots use one canonical
representation for real and complex batch state.

Released recoverable ALPS 3.0.0 linear bin histories migrate offline using
[`alps-hdf5-convert --alea-batches GROUP`](../../tools/hdf5/README.md). Conversion
keeps original bin sums and weights, including partial bins, and recomputes native
analysis uncertainty. Summary-only or nonlinear results cannot recover missing
joint covariance; analysis bins cannot reconstruct an accumulator's merge cursor.

Native `simplemc`/`spinmc` record the blocking hierarchy alongside batch evidence
for every raw measurement, and checkpoint both exactly. Their result files retain
each chain's ordered-series evidence and autocorrelation result independently.
The diagnostic tutorials use these native records directly. Keeping the batch
capacity at least the production length retains the complete individual sample
stream for arbitrary lag analysis or NumPy histograms; the default bounded
capacity retains bin averages, which cannot reconstruct within-bin fluctuations.

## Consolidation coverage and remaining migration

| Capability | Current path |
| --- | --- |
| Real/complex mean, variance, covariance, autocorrelation and batches | Native C++ and Python; circular and elliptic complex strategies |
| Exact accumulator restart, including partial batches/hierarchies | Native codecs kinds 6–10; staged validation |
| Weighted independent-run combination and MPI reduction | Native centered moment algebra; individual chain bins retained |
| Correlated nonlinear estimates and signed ratios | Joint evidence with jackknife or covariance propagation |
| Thermodynamic variances and energy cross-responses | Centered physical moments within spin-engine batches; explicit unavailable reasons |
| Ordered samples, blocking and autocorrelation diagnostics | Separate per-chain native records; full samples when batch capacity permits |
| Histogram analysis | NumPy on retained individual samples; released histograms preserved by conversion |
| Mixed-family C++ joins | Preserve the common available evidence; summary joins assume independence |
| Spin application restart and MPI | Both supported RNGs; extend production; change process count on restart |
| Released statistical files | Offline profile conversion; recoverable histories can become native analysis batches |
| Released scheduler/Parapack application checkpoints | Finish runs with ALPS 3.0; results-only conversion is the supported migration boundary |

The imported `sampling_prop` and `bootstrap_prop` were declarations marked
unimplemented, with no transform implementation. They have been removed from the
API and capability table; this removes no working propagation algorithm.

`mcbase`, CT-INT, CT-HYB, Hirsch-Fye, `simplemc`, `spinmc`, `qwl`,
`dirloop_sse`, `worm`, and `loop` use native estimators. Loop replica ladders
retain separate per-temperature results and checkpoint temperature feedback,
walker permutations, and every physical configuration. See the [quantum solver run guide](../../apps/qmc/README.md)
for checkpoint state, signed correlations and particle-number tuning. The impurity solver analysis files
still lack complete physical solver state and do not provide restart. These
boundaries are explicit; consolidation is not complete until the remaining
clients are ported. Released runs finish with ALPS 3.0; only their results are converted.

QWL uses mean estimators for its coefficient/histogram snapshots and native
batch/autocorrelation estimators for traversal times. It shares the native
runner's serial/MPI chains, RNGs, publication and restart. Exact checkpoints
retain refinement and multicanonical state; completed production can be extended.
Thermodynamic evaluation uses each chain's final corrected coefficients before
averaging curves. See [QWL](../../apps/qmc/qwl/README.md) for window normalization
and the released per-run result converter.

## Removed post-processing APIs

The retired `mcdata.hpp`, `mcanalyze.hpp` and `value_with_error.hpp` headers are
no longer installed or included by `alea.h`. Their application and example
consumers have migrated. The remaining legacy observable producers are still
present until the scheduler/model clients are migrated; this is not completion
of the whole ALEA consolidation.

Use native result means, variances, covariance, errors and `merge` for statistical
analysis. Arithmetic and component selection use `transform`; joint batch
samples retain correlations for nonlinear jackknife propagation. `join` of
summary results assumes independence, whereas aligned batch results retain
joint evidence. Do not reproduce the old practice of treating correlated
quantities as independent error bars.

Raw chronological histories are ordinary arrays, separate from compressed
statistical results. Python's array analysis functions retain lag selection,
exponential fits, integrated times, cuts and running means. The C++ examples in
`tutorials/00-examples/alea` demonstrate the same numerical operations with
STL/Eigen arrays and native moments. They do not introduce another time-series
container. Choose a native accumulator's base batch size before inserting
samples; changing its configuration resets its data. To discard startup samples
or regroup a raw history, slice/group the chronological input and recompute.
A pooled result is not a time series: never assume its batch storage order is
chronological or infer missing samples from bin means.

For published estimates lacking statistical evidence, Python's `ReportedEstimate`
preserves the reported fields without inventing covariance, weights or restart
state; `FloatWithError` supports explicitly independent error arithmetic.
Released archives use the offline `alps-hdf5-convert` profiles (`--alea-batches`,
`--alea-autocorr`, `--alea-summary`). Original fields remain under `legacy/`.
Converted analysis results do not constitute physical simulation checkpoints.

## Remaining consumer boundary

The deterministic DMRG and full/sparse diagonalization applications use generic
scheduler `Task`/`DiagTask` and dispatch, not `MCRun`, `MCSimulation`, or
`ObservableSet`. Their generic scheduler machinery can remain when the Monte
Carlo statistics path is removed. Rewriting that task framework is not required
by ALEA consolidation. The standalone FQHE kernel only needs the HDF5 archive.

The native `loop` executable owns measurement evaluation in `analysis.hpp` and
its current measurement implementation. The former `looper/evaluator.h` and
`evaluator_impl.h` wrappers had no callers, referred to retired evaluator
selectors, and have been removed. No measurement algorithm was removed with them.

Two older source directories require separate decisions before deletion:

- `src/apps/qmc/sse` was already excluded by the released v3.0.0 QMC CMake file
  (`1950cc6f682d7c4c1deae8b816f283857b1819d1`). Its old local CMake file names
  `dirloop_sse_v1` and `dirloop_sse_evaluate`. The built `dirloop_sse` is the
  native `sse4` implementation. Being unbuilt is not a feature-parity proof.
- `src/apps/qmc/sse2` has no CMake target in that release or this tree. Its
  `SIMULATION_PHASE=2` workflow loads external `LOGG_FILENAME` weights, performs
  optimized-ensemble sweeps, and records up-walker histograms for iterative
  weight optimization. The current QWL implementation has Wang–Landau and
  multicanonical sampling but does not implement that external-weight workflow.
  Retain these sources until this unique capability is resolved; do not label
  them a redundant QWL copy or claim complete optional-algorithm parity.

The remaining migration work includes live Parapack/Monte Carlo scheduler
consumers. Released physical checkpoints remain with ALPS 3.0. Removing dead wrappers
or unused includes does not establish completion of those migrations.

### Parapack example migration requirements

The `tutorials/00-examples/parapack` targets are not all redundant copies of
the native solvers. Inspection of their workers and adapters establishes the
following requirements. Except for the completed `ising_single` and
`ising_multiple` ports described below, these remain outstanding work; old
implementations are not correctness oracles.

| Command | Functionality to preserve or establish an equivalent for |
| --- | --- |
| `ising_single` | Native port completed: graph-colored heat-bath updates, deterministic OpenMP within a lattice, extensive moments and centered heat capacity, and resumable warm-start temperature scans. Its shared kernel still supplies the temporary legacy exchange adapter. |
| `ising_multiple` | Native port completed: MPI spatial decomposition of **one** chain, ghost spins and global moments, sharing the serial model's statistics and scans. Odd/uneven ring partitions and exact restart with different rank counts are supported; arbitrary graphs retain a serial fallback. Spatial ranks do not multiply samples. |
| `ising` | Lattice-library heat-bath simulation with bond-type `J`, `J0`, … couplings, normalized magnetization, extensive energy, Binder/heat-capacity analysis, and replica exchange. |
| `heisenberg` | Unit-vector Metropolis updates with bond-type couplings; both vector-magnitude and z-component second/fourth moments and Binder ratios; replica exchange. |
| `loop_single` | Continuous-time quantum loop example with energy, staggered magnetization and uniform/staggered susceptibility estimators. Compare its model, normalization and disorder inputs against native `loop` before consolidating. |
| `exchange` | Classical and quantum workers; serial replica ladders, MPI-distributed replicas, and nested MPI replica/spatial decomposition. Temperature-ladder optimization and exchange diagnostics also belong to this interface. |
| `wanglandau` | Classical Ising **energy** density-of-states learning, fixed-weight microcanonical measurements, and reweighting over temperature, including entropy/reference normalization. The native quantum QWL expansion-order workflow is not an equivalent implementation. |

The spatial Ising worker's acceptance weight and energy sign have been corrected
against the ring Hamiltonian, and its halo transfers use `sendrecv`. It rejects
inconsistent model inputs collectively before communication. Direct worker tests
at one, two and three MPI ranks cover odd and uneven partitions, canonical
ferromagnetic/antiferromagnetic moments, physical bounds, exact restart and
invalid-input consensus, with synchronous sends to expose buffering-dependent
deadlocks. The native command additionally tests partition-independent results,
cross-rank restart of unfinished warmups/statistical bins, both RNGs, temperature
scans, heat capacity and Binder analysis, collective stopping and root-only
publication. Proposals and results also remain identical under OpenMP. Native
spatial preflight checks run settings and input-file bytes across ranks before
physical collectives; constructors and checkpoint loading stay local. The legacy
MPI scheduler remains a separate validation task; its
historical golden files were not scientific references and have been removed.

The authoritative registrations are each directory's `.C` files; the behavior
is in the worker headers and `alps/parapack/{temperature_scan,exchange,exchange_multi,wanglandau}.h`.
The former `temperature_scan_adaptor` retains the worker's physical state
between temperatures, resets measurements after thermalization, and checkpoints
the stage and counters. Independent TOML jobs alone do not preserve this
warm-start workflow. The native `single_ising` implementation retains the spin/RNG state and
statistics per temperature, including separate initial thermalization. Stage
position is derived from its completed-sweep count and immutable stage lengths.
It uses the existing native ensemble archive layout and Python loader; it does
not introduce another scan scheduler or another result format.

Exchange includes explicit temperature/inverse-temperature sets, regular grids,
exchange intervals and randomized ordering, plus rate and population ladder
optimization. Native `<alps/mc/temperature_grid.hpp>` implements both
optimization methods with worker-provided weight laws.
`<alps/mc/replica_exchange.hpp>` now owns the assignment of walkers to temperature
slots, RNG, feedback schedule and restart state, and the loop application uses
this shared implementation. Applications supply physical updates, weights and
statistical recording; signed classical energies are supported without imposing
quantum expansion-order constraints. Acceptance/population diagnostics refer to
temperature slots, while round-trip diagnostics retain physical walker identity.
Do not introduce another exchange engine or restore the legacy observable
framework to reuse the old adapters. The existing native grid requires positive
finite inverse temperatures; the classical beta-zero endpoint needs explicit
assessment during its port, especially for temperature-coordinate feedback.
The generic `alps::mc` runner supplies execution and transport; an application
still has to implement its ensemble's physical state and measurements.

Port and validate these behaviors before deleting their workers or the legacy
Parapack MC framework. Validation must distinguish independent-chain MPI,
spatial MPI, replica exchange, and nested decomposition, and must cover physical
results as well as exact native continuation. Old formulas and golden output
are evidence to inspect, not automatic correctness oracles. The scheduler
Ising migration, for example, exposed an even-sweep parity trap when all flips
were accepted; symmetric proposals that may keep the current spin now avoid
that trap, with a two-site high-temperature regression.
