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
are stored; its other fields are derived. Merging unrelated time series into a
resumable batch accumulator is unsupported; variance and covariance estimators
retain weighted result merging. Other accumulator types have no
checkpoint overload: saving a result is not a resumable checkpoint. Autocorrelation
results retain their reduction API; inserting a result into a live autocorrelation
accumulator is unsupported because it cannot reconstruct the partial hierarchy.

Python exposes `pyalps.alea.BatchAccumulator` and `ComplexBatchAccumulator` with
the same native estimators and checkpoint codecs. Samples are scalar or 1D numeric
arrays; means and errors are always vectors, including one-component results.
`result()` exports a snapshot. `save(archive, path)` / `load(archive, path)` store
and replace complete state; `BatchAccumulator.read(archive, path)` and
`BatchResult.read(archive, path)` construct only after a successful read. Returned
NumPy arrays own their data. Batch sums have `[slots, components]` axes.
`VarianceResult` and `ComplexVarianceResult` read ordinary componentwise
variance results through the same native codec; they expose means, errors,
variance, counts and squared weights without introducing another accumulator API.

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
