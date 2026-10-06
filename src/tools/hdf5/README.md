# HDF5 archive conversion

`alps-hdf5-convert` moves legacy ALPS leaf encodings into the ordinary datatype
mappings used by h5py and HighFive 3.3. It requires Python 3.11+, NumPy and h5py;
it does not import pyalps or load an ALPS library. The SDK installs the script in
its `tools` component. It can also run directly from this checkout:

```sh
python -m pip install numpy h5py
python src/tools/hdf5/convert.py old.h5 converted.h5
# Or, after installing the SDK tools:
alps-hdf5-convert old.h5 converted.h5
```

The source is opened read-only. The destination must not exist. Conversion uses
a temporary file beside the destination and publishes it only on success; a
failure removes the temporary file. There is no in-place or overwrite mode.
The destination filesystem must support hard links for atomic publication.

## Datatype mapping

| Input | Output |
| --- | --- |
| Real/integer components with a trailing dimension of two and scalar `__complex__` marker | HDF5 compound with `r` and `i` members, removing that trailing dimension |
| Signed bytes marked `__alps_type__ = "bool"` | HDF5 enum with `FALSE` and `TRUE` members |
| Signed bytes marked `__alps_type__ = "int8"` | Ordinary signed bytes |
| Explicitly selected ALPS 3.0.0 flat parameter groups | Indexed `alps.params.v2` names and native values |
| Explicitly selected ALPS 3.0.0 ALEA observable/result groups | Native Boolean flags and schema-established empty array extents |
| Explicitly selected ALPS 3.0.0 complete linear bins or equal-weight jackknife histories | Native ALEA batch analysis results; source statistical fields retained as provenance |
| Explicitly selected ALPSCore 2.3.3 ALEA result groups | Canonical primitives and versioned modern ALEA result metadata |
| Other datasets and attributes | Their existing datatypes and values |

Generic complex leaf conversion preserves component width and byte order,
including integer components. The parameter profile instead canonicalizes
numeric widths as described below. Floating complex values are readable as NumPy complex arrays.
Components are assigned separately to preserve infinities, NaNs and signed zero.
Converted Boolean payloads and fill values must contain only zero or one,
including raw codes in selected existing Boolean enums. Consumed ALPS markers
are removed; generic leaf conversion adds no replacement type markers or envelope.
The parameter profile writes its domain checkpoint version and indexed entries.
The same conversions apply to attributes using `__complex__:NAME` and
`__alps_type__:NAME` sibling markers. Malformed, orphaned, unknown and conflicting
markers fail conversion.

Old ALPS and ALPSCore files do not distinguish Boolean from signed-byte storage.
Unmarked signed bytes remain integers and are listed in the conversion report.
Values of zero and one alone do not establish Boolean intent. Declare known
Boolean locations explicitly, including locations reached through hard links:

```sh
alps-hdf5-convert old.h5 converted.h5 \
  --boolean /measurements/enabled \
  --boolean-attribute /simulation converged
```

Both options can be repeated. They cannot override an explicit integer or complex
marker, and misspelled paths fail conversion.

## Released parameter and ALEA profiles

The supported release boundary is files produced by the
[ALPS v3.0.0 writers](https://github.com/ALPSim/ALPS/tree/v3.0.0), revision
`1950cc6f682d7c4c1deae8b816f283857b1819d1`. Select each parameter group and
each individual numeric scalar/vector observable or result explicitly:

```sh
alps-hdf5-convert old.h5 converted.h5 \
  --parameters /parameters \
  --boolean /parameters/ENABLED \
  --alea /simulation/results/Energy \
  --alea /simulation/results/Magnetization
```

Both profiles are repeatable; selections cannot overlap. They never infer a
domain from a group name. Unselected groups retain their scientific schema.

`--parameters` accepts the flat key/dataset groups written by the released
NGS `params` and legacy `Parameters` APIs. It stores each key in
`entries/{n}/name` and its native payload in `entries/{n}/value`, with
`format = "alps.params.v2"`. Datatype and scalar/vector rank determine the value
type; no per-entry type tag is written. Integers are widened to signed/unsigned
64-bit storage, floating values to 64-bit and floating complex values to two
64-bit compound members. Unsupported types or higher ranks are rejected.
Known primitive NULL vectors become ranked empty arrays. Unmarked signed-byte
parameters need an explicit `--boolean` declaration because they could also
come from Python integer storage. String values, including unresolved parameter
expressions, stay text; they are never evaluated or parsed as numbers.

Released parameter writers did not encode keys. Literal `&#47;` and `&#38;`
text therefore remains literal text and is not decoded. An actual `/` became
nested HDF5 paths without recording the original flat-key intent; the profile
rejects nested groups instead of inventing names. Parameter names and values
are validated as UTF-8 without interior NULs. Old string writers placed UTF-8
bytes under an ASCII declaration; selected parameter text is rewritten with a
UTF-8 datatype while retaining supported layout, filters and fill values.
Group and payload attributes and hard-link identity are retained. Soft links
whose target would change after indexing cause conversion to fail.

`--alea` covers the numeric scalar/vector field layouts from `mcdata`,
`SimpleObservableData`, `SimpleBinning` and `DetailedBinning`. It converts their
known scalar `cannotrebin`, `changed` and `nonlinearoperations` attributes to
Boolean enums. These are distinct domain flags; no flag is renamed. Numeric
statistics, convergence flags, bin counts, attributes and stored bin values
retain their meaning. In particular, the profile does not reinterpret sums
as means or infer a missing bin size.

A non-NULL mean, sum, partial bin or bin array establishes the scientific
element shape. Empty scalar bin arrays then have shape `[0]`; empty vector bin
arrays have shape `[0, components]`. The original writer sometimes stored an
empty outer vector as INT NULL, losing its element datatype as well as shape.
Value bins require a non-NULL value-bin/partial-bin exemplar to recover that
type: floating means alone are insufficient because integer measurements can
have floating averages. Logarithmic/jackknife result bins use a result-type
exemplar. Without the required evidence, the profile rejects the file.
Scalar observable labels remain scalar text, and vector labels and logarithmic
counts retain rank one. Generic NULL fields
outside the selected schema remain NULL and are reported.

The compact fixture `tests/cli/fixtures/alps-v3.0.0-profiles.h5` and adjacent
JSON record the pinned writers, source revision and fixture hash. It is an
explicitly labeled h5py reconstruction of those writer contracts, not output
from a compiled release. Tests cover both recoverable empty bins and a genuine
writer layout whose empty vector lost the necessary element evidence.

Private `alps.params.v1` checkpoints created during development of this branch
are not an official release schema and have no automatic upgrade profile.
Other scientific result schemas require an explicit converter profile.

For analysis with native ALEA, use `--alea-batches GROUP` instead of `--alea`:

```sh
alps-hdf5-convert old.h5 converted.h5 \
  --alea-batches /simulation/results/Energy
```

This profile recovers floating real/complex scalar or vector linear bins from
released `DetailedBinning`, untransformed `SimpleObservableData` and `MCData`.
It accounts for their distinct sum/mean conventions and retains the actual
partial-bin weight. Vector bins preserve cross-component covariance. Native
kind-5 means and errors are recalculated from these weighted bins; the uncertainty
can differ from a published legacy estimate.

For transformed `MCData`/`SimpleObservableData` with a complete equal-weight
jackknife history, the same option constructs pseudovalues
`p_i = N*j_0 - (N-1)*j_i`, where `j_0` is the full transformed estimate and
`j_i` are the N leave-one-out estimates. Their native weighted mean and error
reproduce the released bias-corrected mean and jackknife error, and vector
histories retain cross-component covariance. The evaluator writer omitted bin
size; for its transformed results, its `count = N*bin_size` contract supplies
the common weight. Partial or discarded nonlinear histories cannot establish
those weights and are rejected.

Every converted batch result keeps the original statistical leaves (including
reported mean/error, raw variance, tau, bins and jackknife estimates) and domain
flags under `GROUP/legacy`. This is archival provenance accessible with h5py,
not another runtime format reader. Native result operations use only the new
batch state. The original file is also left unchanged. The reserved `legacy`
name must be absent on input; collisions fail rather than overwrite user data.
Summary-only results, missing or inconsistent histories, and aliases to moved
statistical fields are rejected. Unrelated metadata remains intact. Separate
scalar summaries cannot supply missing joint covariance.

For raw `DetailedBinning` histories, the bins and partial bin must reproduce the
recorded total `sum`. The writer added each sample to both, so recursive
summation bounds their difference by `(count + batches) * eps * sum|x|`. The
recorded `sum2` bounds `sum|x|` by `sqrt(count * sum2)`, which also holds when
samples cancel. Without `sum2`, the bin magnitudes are used. A larger
difference means that bins are missing, and the conversion fails.

The output is an analysis result. It cannot supply the merge cursor, RNG or
simulation configuration needed to resume an old NGS checkpoint. New NGS
simulations use native kind-6 batch checkpoints directly and never infer restart
state from legacy results.

## Whole result groups

Released task files store one group per observable under `/simulation/results`.
To convert every observable from the best available evidence, select the
containing group:

```sh
alps-hdf5-convert task.out.h5 native.h5 \
  --parameters /parameters --alea-results /simulation/results
```

Observables with a linear bin history are converted as with `--alea-batches`.
Observables without bins, such as constants that Parapack records without
binning, keep their reported statistics as with `--alea-summary`. The command
prints the choice for each observable. A history that exists but is
inconsistent fails the conversion and never falls back to its summary. Every
child of the selected group must be an observable group.

Task files written by compiled spinmc and loop programs from the 3.0 line test
this path: `tests/cli/fixtures/alps-master-639458499-{spinmc,loop}.h5`. Their
JSON files record the producing revision, parameters and commands.

## Binning diagnostics

`--alea-autocorr GROUP` converts a nonempty real scalar/vector `SimpleBinning`
or `DetailedBinning` hierarchy to the existing native autocorrelation result:

```sh
alps-hdf5-convert old.h5 diagnostics.h5 --alea-autocorr /simulation/results/Energy
```

The released writer stores completed-sample sums, sums of squared **bin means**,
and the number of complete bins at each power-of-two width. The difference
between the total sum and a level's completed sum recovers the unfinished bin.
Conversion therefore retains all samples at every level with their actual bin
weights. Native unbiased weighted variances and errors are recomputed, including
single-bin infinite uncertainty. Every stored level remains available; the old
plotting reader's unconditional removal of the last seven levels is gone.
Native error selection uses its usual 1024-observation threshold, so published
legacy errors and convergence flags may differ. They remain archived under
`legacy/`, alongside all original statistical fields.

`loadBinningAnalysis` reads the converted native levels without decoding legacy
moments at runtime. Missing/inconsistent counts, nonfinite moments, transformed
histories and reserved-name collisions fail conversion. Complex squared moments
cannot supply the circular variances required here and are rejected. Tiny negative
centered moments within floating-point roundoff are clamped to zero; larger
inconsistencies fail. Conversion cannot recover precision already lost when the
old writer accumulated raw moments. This creates analysis diagnostics, not an
accumulator cursor or simulation restart. Select `--alea-batches` instead when
the desired output is a joint batch result; profile selections cannot overlap.

## Reported estimates without recoverable histories

Use `--alea-summary GROUP` when the recorded statistics must be retained without
interpreting incomplete bins, or when only a published estimate remains:

```sh
alps-hdf5-convert old.h5 reported.h5 --alea-summary /simulation/results/Energy
```

The output has `format="alps.reported-estimate.v1"`, a uint64 `count`, and
one-dimensional component arrays at `mean/value` and, when supplied, `mean/error`,
`variance/value`, `tau/value` and `mean/error_convergence`. Floating statistics
use float64 (complex128 means are supported); convergence flags use int64.
Optional fields remain absent rather than becoming zero. Undefined numerical
estimates remain as reported. Original statistical fields and histories remain
under `legacy/`, using the same provenance preservation as `--alea-batches`.
There is no native estimator `kind`, inferred covariance, effective sample count,
reconstructed batch or resumable state.

Python `alea.read_result` and `loadMeasurements` read this explicit format as a
`ReportedEstimate`; `saveMeasurements` can save its reported fields again.
The mean and variance commands also accept it, but variance requires an actual
reported variance. A published standard error alone cannot determine it.
Native pooling rejects these records: they lack the sampling evidence required
by that algorithm. Use recoverable native batches when subsequent correlated
analysis or native pooling is needed. When resaving reported estimates, the
record's fields are saved; archive-level provenance stays in the converted file.

## Released ALPSCore ALEA results

For modern ALEA results written by ALPSCore 2.3.3, select the estimator explicitly:

```sh
alps-hdf5-convert core.h5 converted.h5 \
  --core-alea covariance /results/EnergyMagnetization \
  --core-alea batch /results/Correlations
```

Kinds are `mean`, `variance`, `covariance`, `autocorr` and `batch`. The profile
validates their field shapes and adds modern ALEA's `version=1` and estimator
`kind` attributes, including variance results nested in autocorrelation levels.
It preserves counts, squared weights, batch sums and per-batch counts, full
covariance and Eigen's physical `[columns, rows]` axes. Complex circular
covariance uses complex compounds; elliptic covariance keeps its real 2×2
operator axes. Results remain results: conversion cannot recover an accumulator's
missing merge cursor or unfinished sampling state. It cannot reconstruct
cross-covariance from separately saved scalar summaries.

`tests/cli/fixtures/alpscore-v2.3.3-alea.h5` was produced by compiled Core writers
from the pinned reference; its adjacent JSON records the producer, fixture hash
and comparison with release `f2ccddc5343bdc2297727f10a1a12d9c717cd0ab`. The
released statistical serializers are unchanged; `batch.cpp` has only an added
`<algorithm>` include. The native migration reader compares all converted
estimators with the same known sample streams, including both complex covariance
conventions.

## Released QWL final results

```sh
alps-hdf5-convert old.out.run1.h5 native.h5 --qwl-sites 40
qwl_evaluate native.h5
```

This whole-file profile converts ALPS 3.0.0 QWL coefficient/histogram means to
native mean results and indexes parameters. Supply the lattice site count; the
profile checks it against the order-zero normalization. Completed per-run files
and single-run summaries are supported. Each final estimate must have count one;
pooled logs cannot recover separate chains, so use individual `.out.runN.h5`
files for multi-run work. Nonzero windows from the released writer are rejected.
Coefficients produced without combinatorial factors receive the corresponding
factorial correction. Other scientific fields, including timing bins and original
sums, remain available. No solver restart state is inferred.

The profile accepts unaliased hard-linked writer layouts and cannot be combined
with other domain selections. The native evaluator consumes the ordinary native
schema; it contains no legacy file reader. The fixture
`tests/cli/fixtures/alps-v3.0.0-qwl.h5` contains actual pre-migration branch writer
output; its metadata records the producer and the comparison with released QWL
measurement writers. It is not represented as output from a compiled v3.0.0 build.

## Application checkpoint conversion boundary

A released result file is not necessarily a complete application checkpoint.
For example, the ALPS 3.0.0 `spinmc` HDF5-enabled writer splits a restart across
its `.out.runN.h5` archive and companion `.out.runN` XDR file:

| State | Released writer and location |
| --- | --- |
| Parameters and RNG name/state | `Worker::save(hdf5::archive&)`: `/parameters`, `/rng` and `/rng/@name` |
| Measurements | `MCRun::save(hdf5::archive&)`: `/simulation/realizations/0/clones/0/results` |
| Total updates and fractional cluster thermalization | `AbstractSpinSim::save(ODump&)`: XDR payload |
| Physical spins | `SpinSim::save(ODump&)`: XDR payload, with model-specific moment representation |

This follows the pinned release sources in `src/alps/scheduler/worker.C`,
`src/alps/scheduler/montecarlo.C` and
`applications/mc/spins/{abstractspinsim,spinsim,ising,potts,on}.h`.
Despite the comment in the released `config.h.in` saying checkpoint data is
HDF5-only, application `ODump` hooks still write physical state. The worker
stream starts with the run marker, reserved integer and format version. Version
400 moves framework state to HDF5; version 310 additionally embeds parameters,
RNG, task information and measurements in XDR. These layouts must not be mixed.

The `--spinmc-state` profile recovers the version-400 physical payload into
ordinary HDF5 datasets while copying the companion HDF5 archive:

```sh
alps-hdf5-convert task.out.run1.h5 recovered.h5 \
  --spinmc-state task.out.run1 --parameters /parameters
```

`/migration/spinmc` contains the model, spins (site × component), total-update
counter, fractional thermalization counter, thermalization sweep counter, and
exact source XDR bytes. Ising Boolean states become ±1; Potts states retain their
color indices; XY, Heisenberg and O(4) retain their vector components. The reader
checks the header, exact payload length and physical spin constraints before
publication. Potts requires a numeric `q` of 3, 4 or 10; unresolved expressions
are rejected. Version 310 is not yet supported. RNG and measurement data remain
in their original archive locations; statistical profiles can be selected in the
same invocation. The source files remain untouched.

The released pair has no shared identifier: callers must supply companions
from the same saved run. Presence and payload checks cannot authenticate that
pairing. The fixtures are compiled C++ protocol reconstructions using OSIRIS,
with source references and hashes in `tests/cli/fixtures/spinmc-state.json`;
they are not represented as complete released application checkpoints.

This profile is state recovery, not native solver restart conversion.
No existing profile translates these physical checkpoints into native solver
checkpoints. Keep both companion files. A complete offline translator must
validate the pair and model representation, carry physical state and RNG into
the native engine, and explicitly preserve the available statistical history.
Missing histories or changed estimator definitions must not become invented
native batches or be silently discarded. A warm start with reset measurements
would be a different operation, not complete checkpoint conversion.

## Older explicit container schema migration

Old pair and numerical matrix groups require an explicit selection; their child
names alone do not establish container intent:

```sh
alps-hdf5-convert old.h5 converted.h5 \
  --pair /simulation/pair \
  --matrix /simulation/matrix
```

Both options are repeatable. A selected pair must contain exactly `first` and
`second`; those links become `0` and `1`. Its group/child attributes, hard-link
aliases and cycles are retained.

A selected matrix must contain exactly `size1`, `size2`, `reserved_size1` and
`values`. Its integer scalar dimensions and flat storage length are checked.
Padded column storage becomes a dataset with shape `[columns, rows]`; reserved
row padding is discarded. Converted complex values retain their compound type.
Explicit dimensions recover empty NULL storage. Group and value attributes are
retained, compression/checksums remain where supported, and chunk dimensions
are recalculated for the matrix shape. Data is copied in bounded tiles.

Matrix groups/fields with hard-link aliases, annotated dimension fields,
conflicting group/value attribute names, extensible storage and unsupported
filters are rejected. Selections must follow hard links and cannot overlap.
Schema migration also rejects changes to soft-link targets or to whether a soft
link resolves. These cases require an explicit mapping of the affected metadata
or links instead of silently discarding them. Without a selection, old pair and
matrix groups remain unchanged. The genuine old matrix fixture is retained at
`tests/cli/fixtures/legacy-alps-matrix.h5`.

## Scope and limits

Group paths, encoded names, unrelated attributes, hard-link identity, cycles and
soft links are preserved. Indexed groups used for ragged containers stay groups:
their original container kind cannot be inferred reliably. Unchanged datasets
are copied by HDF5; transformed datasets are processed in bounded tiles.
Chunking, extendible dimensions, gzip/LZF compression, shuffle and checksums are
retained where compatible. Complex scalar conversion produces a scalar dataset,
which cannot retain chunking or compression.

Generic NULL dataspaces remain NULL and are reported. Older empty-container
writers lost the original array rank; the converter does not invent one.
HighFive 3.3's ordinary vector read does not accept NULL as an empty vector.
Ranked zero-length arrays retain their shape and can be read normally.

Files with user blocks, external links/storage, virtual datasets, named datatypes,
object/region references (including dimension-scale references), and unsupported
filters on transformed datasets are rejected. These cases require a deliberate
mapping; copying references or moving relative external paths could change data.
Only a single, self-contained input file is supported.

Other scientific checkpoint schemas and application version values remain
unchanged. The tool does not recover original types from textual parameters,
reinterpret random-generator state, merge statistics, or translate solver versions.

Legacy conversion belongs here. New runtime code should use the selected
HighFive/HDF5 mappings directly; adding an automatic old-format reader there
would reintroduce the compatibility burden this tool is intended to remove.

## Scientific capability boundary

Conversion and estimator reconstruction are different operations. The following
inventory is the gate for retiring the remaining legacy measurement clients:

| Released evidence | Conversion / analysis capability | Information that cannot be inferred |
| --- | --- | --- |
| Mean/error/variance summaries and convergence flags | `--alea` preserves values, widths and their distinct meanings, including nonlinear results | Joint covariance, effective independent sample count, or jackknife bins |
| Histogram counts with `min`, `max`, `stepsize` | Generic conversion preserves this ordinary HDF5 layout; NumPy can use the counts directly | Original sample order or values within a histogram cell |
| Linear histories including a partial bin | `--alea` retains the history; `--alea-batches` produces native weighted analysis evidence | Individual observations inside a bin |
| Logarithmic histories | `--alea` retains all stored logarithmic sums, squared sums and counts | Missing levels, discarded samples, or a missing unfinished bin |
| Core 2.3.3 native result families | `--core-alea` retains all saved statistical evidence | Unsaved accumulator state |
| Released application checkpoint | Generic conversion preserves unselected state fields; it does **not** translate the application's restart protocol | A new solver's configuration/update state from an analysis summary |

The histogram contract is `HistogramObservable<T>::save` in ALPS v3.0.0
`src/alps/alea/histogram.h`, at the release revision above: unsigned 32-bit cell
counts, unsigned 64-bit total count, and range/step attributes. No histogram
format adapter or new runtime class is needed. Tests independently reconstruct
this contract and verify preservation alongside summaries and histories.

Keep original released checkpoints when continuation with the released program
is needed. The new spin engines checkpoint their complete native state, and can
continue those files across serial/MPI process counts; a released scheduler or
Parapack checkpoint is not yet convertible into that new application state.
This remains a migration gap, and blocks a claim of complete checkpoint-format
coverage. Do not remove an old checkpoint reader merely because statistical
result conversion passes. CT-INT, CT-HYB and Hirsch-Fye analysis files omit solver
configurations; a converter cannot turn them into resumable solver checkpoints.
