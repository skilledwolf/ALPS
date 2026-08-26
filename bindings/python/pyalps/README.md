# pyalps

Python applications and libraries for the Algorithms and Libraries for
Physics Simulations (ALPS) project. Binary wheels are available from PyPI:

```sh
python -m pip install pyalps
```

Install `pyalps[plot]` to use the Matplotlib plotting helpers.
Install `pyalps[mpi]` for the mpi4py-backed `pyalps.mpi` compatibility layer.

The bindings are built as a standalone `scikit-build-core` project using
nanobind. A source build requires Python 3.10 or newer, CMake 3.22 or newer,
Ninja, a C++17 compiler, BLAS/LAPACK, HDF5, and an installed ALPS C++ SDK.
Point `ALPS_DIR` at the SDK's `share/alps` package directory.

The `wheel-deps` CMake preset builds the SDK exactly as the wheel CI does.
From the repository root:

```sh
cmake --preset wheel-deps
cmake --build --preset wheel-deps

ALPS_DIR="$PWD/_build/wheel-deps/install/share/alps" \
  python -m build --wheel bindings/python/pyalps
```

The wheel is written to `bindings/python/pyalps/dist` and can be installed
with `python -m pip install`. With ccache installed, configure with
`cmake --preset wheel-deps -DCMAKE_CXX_COMPILER_LAUNCHER=ccache` and set
`CMAKE_ARGS="-DCMAKE_CXX_COMPILER_LAUNCHER=ccache"` for the wheel build to
speed up rebuilds.

`PYALPS_BUILD_APPLICATIONS=ON` is the default and preserves the MaxEnt,
DWA, CT-HYB, and CT-INT extension modules. Set it to `OFF` through CMake
configuration for a smaller core-only developer build.

## Free-threading and stable-ABI policy

pyalps ships **one stable-ABI (`cp310-abi3`) wheel per platform** that
covers CPython 3.10 and newer. The extension modules build in nanobind's
split mode (`BACKEND_MODULE nanobind_backend`): they contain only the
tiny nanobind frontend, compile under `Py_LIMITED_API` (3.10 floor), and
resolve the compiled nanobind runtime at import time from the
`nanobind-backend` package, which is a runtime dependency of pyalps.
Consequences:

- **Binding code must stay limited-API clean.** In particular,
  `PyTypeObject` is opaque — type-name dispatch goes through
  `alps::python::type_fullname()` (cpp/numpy_compat.hpp) instead of
  `tp_name`. Violations fail at compile time, so a successful CI build
  is the enforcement.
- **Downstream extensions** that need bound-type identity with pyalps
  (e.g. `mcbase` subclasses) must also build in split mode against the
  same backend module — see `tutorials/ngs/5_export_python`. Extensions
  only share nanobind type bindings when they share a backend.
- **Free-threading (3.13t/3.14t):** still deliberately unsupported; abi3
  wheels do not install on free-threaded interpreters. The ALPS C++
  library relies on the GIL as its lock around shared state
  (`mcobservable`'s reference-count table, the `alps::ngs::signal`
  singleton, `mcdata`'s lazily-computed statistics). Do not add
  `FREE_THREADED` to `nanobind_add_module` without first making that
  state thread-safe.
