<!--
Copyright (c) 2026 Tobias Wolf
SPDX-License-Identifier: MIT

Adapted from skilledwolf/alpsx (private repository), branch
agent/v3-execution-and-api-audit at d37e51227,
python/qdk-chemistry-alps/docs/install.md. Changes: rewritten for the ALPS SDK
and pyalps of this repository (CMake SDK build, pyalps bindings, no Pixi
environment, no NumPy override); the RdmContainer build steps are kept.
-->

# Installation

The package needs three things in one Python environment: the ALPS SDK, pyalps built against it with its `models`, `ed` and `classic_dmrg` bindings, and qdk-chemistry. pyalps is not a declared pip dependency, because the bindings come from a source build against the SDK of the same checkout. The `pyalps` wheel on PyPI predates these bindings.

## Prerequisites

- Python 3.12 or newer (pyalps itself supports 3.10 and newer).
- The ALPS SDK prerequisites: CMake 3.27+, Ninja, a C++17 compiler, Boost 1.76+, LP64 BLAS/LAPACK and HDF5 1.10.5+. The [contributor setup](../../../CONTRIBUTING.md#build) lists the supported stacks.
- qdk-chemistry 2.1 or newer. Released 2.2.1 from PyPI works; only the RDM-only returned wavefunction needs [a build providing `RdmContainer`](#a-qdk-chemistry-that-provides-rdmcontainer).
- Optional: `pyscf` for the molecular examples and reference tests, `matplotlib` for the figures (`--output` on the examples, `tools/generate_figures.py`).

## 1. Build and install the SDK

From the ALPS repository root:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PWD/_alps-install" \
  -DALPS_ENABLE_MPI=OFF
cmake --build build --parallel
cmake --install build
```

`ALPS_ENABLE_CLASSIC_DMRG` is on by default; `alps_classic_dmrg` needs it. The applications are not needed by this package. Use a new build directory when changing compilers.

## 2. Install pyalps against the SDK

```bash
python3.13 -m venv .venv
source .venv/bin/activate
ALPS_DIR="$PWD/_alps-install/share/alps" PYALPS_BUNDLE_APPLICATIONS=OFF \
  python -m pip install ./python/pyalps
```

`PYALPS_BUNDLE_APPLICATIONS=OFF` builds the bindings only, against the installed SDK; `PYALPS_BUILD_MODELS` (on by default) builds the `models`, `ed` and `classic_dmrg` bindings. The [pyalps README](../../pyalps/README.md#building-from-source) covers wheels, editable installs and the build options.

## 3. Install qdk-chemistry and this package

```bash
python -m pip install 'qdk-chemistry>=2.1' pyscf pytest
python -m pip install --no-deps -e python/qdk-chemistry-alps
```

Without pyalps the import of `qdk_chemistry_alps` fails with this guidance.

## 4. Verify

```bash
python -c "
import qdk_chemistry_alps
from qdk_chemistry.algorithms import available
print(sorted(n for n in available('multi_configuration_calculator') if n.startswith('alps')))
"
python python/qdk-chemistry-alps/examples/quickstart.py
python -m pytest python/qdk-chemistry-alps --require-solver-stack -ra
```

Expect `['alps_classic_dmrg', 'alps_ed']`, and ED and DMRG energies of `-3.0925653195` that agree within `1e-8`. `--require-solver-stack` fails instead of skipping when pyalps, qdk-chemistry or the classic DMRG bindings are missing, so a packaging-only run cannot pass for a numerical one. Inspect the skip report: tests marked `pyscf`, `local_qdk` and `qsharp` skip without their prerequisites.

| Capability | Needed for | Test marker |
| --- | --- | --- |
| `pyscf` | Molecular examples and reference energies | `pyscf` |
| A qdk-chemistry providing `RdmContainer` | RDM-only wavefunctions past the CI-vector limit | `local_qdk` |
| pyalps `classic_dmrg` bindings | `alps_classic_dmrg` | `classic_dmrg` |
| The `qdk` package's Q# simulator, which qdk-chemistry installs | The two quench examples | `qsharp` |

Without `RdmContainer`, the DMRG calculator still exposes the RDMs through `make_rdms()` and returns state vectors, but raises `RuntimeError` naming the class before solving when its settings call for an `RdmContainer`. See [Work past the CI limit](how-to/large-active-spaces.md).

## A qdk-chemistry that provides `RdmContainer`

Released qdk-chemistry 2.2.1 does not provide this container. Use a branch that implements it; a newer version number alone does not imply support. From such a checkout, in the environment above:

```bash
cmake -S cpp -B cpp/build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PWD/_install"
cmake --build cpp/build --target install --parallel
CMAKE_PREFIX_PATH="$PWD/_install" python -m pip install --no-deps ./python
```

Setting the install prefix at configure time also directs bundled dependencies to that prefix.

## Troubleshooting

**Missing native function or attribute** in `pyalps.models`, `pyalps.ed` or `pyalps.classic_dmrg`: the Python sources and the compiled bindings come from different revisions. Reinstall the SDK and pyalps from the same checkout and restart Python.

**CMake cannot find ALPS:** point `ALPS_DIR` at the installed SDK's `share/alps` directory. A CMake build tree is not an installed SDK.

**Imports resolve to another pyalps:** check `python -c "import pyalps; print(pyalps.__file__)"`, uninstall the stray copy in that environment and reinstall. Do not work around a binding mismatch with `PYTHONPATH`.
