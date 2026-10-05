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
| Explicitly selected ALPS 3.0.0 complete linear bin histories | Native ALEA batch analysis results with original sums/counts and recomputed uncertainty |
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
can differ from a published legacy estimate. Summary-only, nonlinear, discarded,
missing or inconsistent histories are rejected, as are aliases to replaced
statistical fields. Unrelated metadata remains intact. Separate scalar summaries
cannot supply missing joint covariance.

The output is an analysis result. It cannot supply the merge cursor, RNG or
simulation configuration needed to resume an old NGS checkpoint. New NGS
simulations use native kind-6 batch checkpoints directly and never infer restart
state from legacy results.

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
