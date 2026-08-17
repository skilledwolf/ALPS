# CI build matrix

`linux.json` and `macos.json` are the single source of truth for the ALPS build
matrix. They are consumed by the reusable workflow
[`build-alps.yml`](../workflows/build-alps.yml), which filters entries by tier:

- **`"tier": "quick"`** — runs on every pull request and push to master
  ([`ci.yml`](../workflows/ci.yml)). These entries cover the *boundaries* of each
  axis: oldest/newest GCC, oldest/newest Clang, newest C++ standard, a combined
  all-minimums job (oldest Boost + oldest Python + C++11), and one macOS
  arm64 + one macOS x86_64 job.
- **`"tier": "full"`** — additionally runs in the weekly sweep and on manual
  dispatch ([`full-matrix.yml`](../workflows/full-matrix.yml)): the interior of
  each version axis (intermediate compilers, Boost, Python, C++ standards) plus
  extra macOS targets. The full tier is a superset — it includes the quick
  entries too.

Rationale: intermediate versions almost never break independently — if the
oldest and newest compiler both pass, the ones in between failing alone is
rare. Boundaries gate PRs; the exhaustive sweep catches ecosystem drift
(new compiler/Boost/Python releases) on a schedule instead of taxing every PR.

Entry fields: `label` (job name + ccache key — keep unique), `tier`, `os`,
`comp_pack` (apt/brew package spec), `c_compiler`/`cxx_compiler`, `c_version`
(Linux only), `cxx_stdlib` (macOS only), `py_version`, `boost_version`
(e.g. `91` = Boost 1.91.0), optional `cxx_standard` (defaults to 14).

Axis sweeps hold everything else at the baseline: ubuntu-24.04, gcc-14,
Python 3.14, Boost 1.91, C++14.
