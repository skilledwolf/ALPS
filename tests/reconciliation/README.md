# ALPS / ALPSCore contract probes

This standalone project compiles the same source against **one installed SDK at
a time**. The runner exchanges temporary HDF5 files between fresh processes.
Never combine the two SDK include paths or runtime implementations in one binary.
The probes characterize differences; a successful comparison does not establish
drop-in compatibility.

The runner covers 20 archive payload cases, params access/conversion behavior,
full params versus dictionary checkpoint loading, and wider scalar checkpoints.
Every typed archive read checks the value. Common and extended params self-checks
verify exact values for every supplied provider; cross-provider differences are
recorded even when loading throws or silently changes a value. `--expect`
additionally rejects any measurement change
from an earlier report. A known defect recorded in a baseline is not a desired
contract: review and update that expectation when fixing it.

ALPS self-reads also require the `alps.params.v2` format. ALPS writes native
compound complex values, Boolean enums and ranked zero-length arrays; the pinned
Core writes older physical encodings. Consequently direct cross-provider reads
fail by design, while every same-provider archive value must round-trip.
ALPSCore's INI/argv/default bookkeeping is not part of the ALPS params layer.

With `--converter`, the runner also migrates each Core params checkpoint with
`alps-hdf5-convert --core-parameters /parameters` and requires ALPS to read the
converted file exactly as its own checkpoint, including the `alps.params.v2`
format, Boolean, unsigned, wide-integer and single-precision values. Archive
payloads other than parameters are converted with the profiles described in the
[converter guide](../../src/tools/hdf5/README.md).

The conversion observations include Boolean-to-integer access, integer-to-float
precision loss at 2^24 + 1, and integer-vector-to-real-vector access. These are
intentional differences: ALPS rejects the first two conversions and converts
homogeneous vectors element by element; the pinned Core permits the first two
conversions and requires matching vector storage types. Observations
characterize each provider and are not assertions that their policies must match.
Core's context-complex and overflowing-string defects remain characterizations,
not desired behavior. ALEA and MPI are outside this comparison.

## Linux x86_64 baseline

[`baseline-linux-x86_64.json`](baseline-linux-x86_64.json) records the comparison
of the ALPS SDK and unpatched Core `7146b9e1f017938a94e5dae35d88467cc5ba7969`,
with the converted-Core check. Both SDKs used GCC 13.3.0, Boost 1.83.0 and HDF5
2.2.0, C++17, Release, shared libraries and MPI disabled; Core also used Eigen
3.4.0 headers with its `accumulators`, `mc`, `gf` and `alea` modules disabled. No
reference source was patched. Earlier macOS baselines measured the
pre-`alps.params.v2` formats and remain in the repository history.

Use a baseline only for its platform and measured contracts. On platforms with
a narrower `long`, the wide-integer fixture uses 42 instead of 2^40, and integer
overflow results depend on the implementation. Run without `--expect` to
characterize another platform or revision, then review the report. The JSON
records measurements, platform and executable hashes; its sibling `.log`
retains exception diagnostics, including why checkpoint loads failed.
Use `--alps-revision` and `--alpscore-revision` to supply known SDK source
revision labels. Omitted labels are recorded as `unknown`; the runner does not
infer a library's revision from the current source checkout or the probe
executable hash. Record dependency versions alongside any new baseline.

## Build and run

Use the already installed ALPS SDK and its compiler/dependencies. When the Core
reference is unavailable, build only the ALPS probe and omit `--alpscore`:

```sh
python tests/reconciliation/compare.py \
  --alps /absolute/path/to/alps-probe/reconciliation_probe \
  --output /absolute/path/to/alps-self-check.json
```

This checks ALPS archive values and common/extended params self-reads. The report
records `same-provider` scope and does not establish ALPSCore interchange.
Native params tests cover the broader canonical storage families.

For a cross-provider run, install only ALPSCore utilities/HDF5/params from an
unmodified checkout of `7146b9e1f017938a94e5dae35d88467cc5ba7969` into a
separate prefix:

```sh
# Set these to existing absolute paths. Keep this scratch tree outside the SDK.
alps_repo="$PWD"
alps_sdk="$alps_repo/_build/default/install"
alps_hdf5="$alps_repo/_build/hdf5-install"
core_source=/absolute/path/to/pinned/ALPSCore
scratch=/absolute/path/to/reconciliation-build

cmake -S "$core_source" -B "$scratch/core-sdk" -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DALPS_CXX_STD=c++17 \
  -DCMAKE_PREFIX_PATH="$alps_hdf5" -DHDF5_ROOT="$alps_hdf5" \
  -DEIGEN3_INCLUDE_DIR=/usr/include/eigen3 -DALPS_INSTALL_EIGEN=OFF \
  -DTesting=OFF -DENABLE_MPI=OFF \
  '-DALPS_MODULES_DISABLE=accumulators;mc;gf;alea' \
  -DCMAKE_INSTALL_PREFIX="$scratch/core-sdk/install"
cmake --build "$scratch/core-sdk" --parallel 2
cmake --install "$scratch/core-sdk"

cmake -S "$alps_repo/tests/reconciliation" -B "$scratch/alps-probe" -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DRECONCILIATION_PROVIDER=ALPS \
  -DALPS_DIR="$alps_sdk/share/alps" -DHDF5_ROOT="$alps_hdf5"
cmake --build "$scratch/alps-probe"

cmake -S "$alps_repo/tests/reconciliation" -B "$scratch/core-probe" -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DRECONCILIATION_PROVIDER=ALPSCore \
  -DCMAKE_PREFIX_PATH="$scratch/core-sdk/install;$alps_hdf5"
cmake --build "$scratch/core-probe"

python "$alps_repo/tests/reconciliation/compare.py" \
  --alps "$scratch/alps-probe/reconciliation_probe" \
  --alpscore "$scratch/core-probe/reconciliation_probe" \
  --converter "$alps_repo/src/tools/hdf5/convert.py" \
  --alpscore-revision 7146b9e1f017938a94e5dae35d88467cc5ba7969 \
  --output "$scratch/results.json"
```

The converter needs h5py and NumPy in the interpreter that runs `compare.py`.
Adjust paths, compiler and executable suffix for other platforms. Pin `HDF5_ROOT`
as well as `CMAKE_PREFIX_PATH` when multiple providers are installed.

The runner removes its generated archive fixtures automatically. Preserve useful
incremental SDKs; remove only the dedicated scratch builds when finished. The
foreign SDK is optional and is not downloaded or built by the main ALPS test suite.

The runner's failure-detection tests require no SDK build:

```sh
python -m pytest tests/reconciliation/test_compare.py -q
```

## Gaps

This is a small contract baseline, not a full format/ABI audit. It does not yet
exercise schema-defined defaults/provenance, scientific-container adapters,
multidimensional shape inspection, compression/family files, concurrent access,
signal/crash cleanup, leak instrumentation, Python providers, MPI, or MaxEnt
scientific results.
