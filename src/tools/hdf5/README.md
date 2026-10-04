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
| Typed `alps.params.v1` dictionaries | Validated `alps.params.v2` dictionaries with native datatypes and ranks |
| Other datasets and attributes | Their existing datatypes and values |

Complex component width and byte order are preserved, including integer
components. Floating complex values are readable as NumPy complex arrays.
Components are assigned separately to preserve infinities, NaNs and signed zero.
Boolean payloads and fill values must contain only zero or one, including raw
codes in existing Boolean enums. Consumed ALPS markers are removed;
the converter adds no replacement type markers or private format envelope.
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

## Explicit container schema migration

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

Explicit `alps.params.v1` dictionaries are upgraded to `alps.params.v2`.
The indexed names and logical types are validated before conversion. Declared
Boolean values resolve otherwise ambiguous signed bytes, and declared vector
types recover NULL values as one-dimensional empty arrays. Scalars must have
scalar rank; vectors must have rank one after complex conversion. Unknown types,
mismatched physical datatypes/ranks, duplicate names and sparse entry indices
are rejected. Names containing `/` or TOML punctuation remain indexed names.
String names and values are validated as UTF-8 without interior NULs. Old writers
stored UTF-8 bytes under an ASCII charset declaration; typed dictionary text is
rewritten with a UTF-8 HDF5 string datatype while retaining supported layout,
compression and fill values.

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
