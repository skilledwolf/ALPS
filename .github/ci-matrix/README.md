# CI build matrix

`linux.json` and `macos.json` are the single source of truth for source-build
coverage. The reusable workflow combines both files and filters entries by tier.

- **`quick`** is the pull-request gate: one all-minimums Linux build, current
  GCC, current Clang with the newest C++ mode, and current macOS arm64.
- **`full`** runs weekly and on demand. It adds the supported Clang floor,
  next compilers, Linux arm64, macOS Intel/current-OS canaries, sanitizers,
  next Python, and a real build without legacy permissive mode.

The quick tier contains four builds; the full tier contains fourteen. We test
support boundaries and meaningful platform combinations rather than every
interior compiler, Boost, Python, and C++ version. Interior versions share the
same language and ABI boundaries and previously consumed most runner capacity
without independently useful signal.

Every Boost entry includes the SHA-256 published beside the archive on
`archives.boost.io`. The build passes the extracted directory as
`Boost_SRC_DIR`, which is the variable consumed by `FindBoostForALPS.cmake`.
The all-minimums combination constrains NumPy below 2 because Boost releases
before 1.87 do not support the NumPy 2 C API.

Only quick entries set `remote_cache: true`. Their ccache archives are restored
on pull requests but saved only from `master`, with one immutable generation per
ISO week. Scheduled compatibility entries intentionally build cold, preventing
the full matrix from creating dozens of disposable caches per run.

Canaries are strict. A scheduled sanitizer, next-toolchain, next-Python, or
no-permissive failure fails the compatibility workflow so it cannot remain a
silent yellow result.
