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
checkpoint overload: saving a result is not a resumable checkpoint.
