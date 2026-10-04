# Modern ALEA statistical core

`ALPS::statistics` contains the Eigen-based ALEA statistical core imported from
ALPSCore revision `7146b9e1f017938a94e5dae35d88467cc5ba7969`. The source and
small `common::ndview` / serialization support headers retain their original
ALPS Collaboration copyright notices and MIT licensing. Core's MPI and stream
codec plugins and package build system are not imported. Existing legacy
ALEA APIs remain with `ALPS::alps` during the measurement pilot.

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

Independent-run result reduction retains every batch and its weight, including
unfinished batches. Runs with different slot counts use disjoint blocks padded
with empty slots; their bins are never summed together. Autocorrelation reduction
retains only levels present in every run, so every retained level contains every
run's samples. Failed reductions preserve the original result. Reduction does
not create a resumable combined time series.

The optional `<alps/alea/mpi.hpp>` reducer takes a borrowed `MPI_Comm`; clients
link `MPI::MPI_CXX` alongside `ALPS::statistics`. The statistics library itself
remains MPI-free. Construction and reductions are collective and must occur in
the same order on every rank. Only the chosen root retains the combined result.
Sample counts use unsigned 64-bit MPI arithmetic.
Custom reducers must implement `reduce(view<uint64_t>)`; rebuild downstream
binaries after this interface change. The combined sample count must fit in
`uint64_t`.

`pyalps.hdf5.save_checkpoint(filename, callback)` calls the existing native
publication helper. The callback receives a `NativeArchive`; retained callback
views close before publication, and a failed save preserves the previous file.
The pure Python Ising tutorial uses this path and checks complete restart against
an uninterrupted spin stream. Both Python and C++ pilots use one canonical
representation for real and complex batch state.
