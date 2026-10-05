# Modern ALEA statistical core

`ALPS::statistics` contains the Eigen-based ALEA statistical core imported from
ALPSCore revision `7146b9e1f017938a94e5dae35d88467cc5ba7969`. The source and
small `common::ndview` / serialization support headers retain their original
ALPS Collaboration copyright notices and MIT licensing. Core's MPI and stream
codec plugins and package build system are not imported. Existing legacy
ALEA APIs remain with `ALPS::alps` while active clients migrate.

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

The NGS `mcbase` framework stores a map of shared real batch accumulators with
explicit component dimensions. Python assigns `BatchAccumulator` handles directly
to `sim.measurements[name]`; a retrieved handle remains valid after replacement
or erasure. Collected results are owning maps/dictionaries of native batch results.
Base checkpoint loading stages parameters, measurements and RNG and validates
all already registered names and dimensions before replacement. Subclasses
override archive-reference hooks for their application state.


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
`mcmpiadapter` checks matching result requests before reducing every entry,
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
