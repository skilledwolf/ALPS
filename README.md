[![Build](https://github.com/ALPSim/ALPS/actions/workflows/build.yml/badge.svg)](https://github.com/ALPSim/ALPS/actions/workflows/build.yml) [![Python wheels](https://github.com/ALPSim/ALPS/actions/workflows/build_wheels.yml/badge.svg)](https://github.com/ALPSim/ALPS/actions/workflows/build_wheels.yml) [![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE.txt)

# ALPS — Algorithms and Libraries for Physics Simulations

The ALPS software package aims to provide a set of well tested, robust, and standardized components for numerical simulations of condensed matter systems, including bosonic, fermionic, and spin systems. They consist of a set of components that are used in state-of-the-art high performance codes.

**Project website:** [alps.comp-phys.org](https://alps.comp-phys.org/)

## Learn ALPS

Start with the [tutorial guide](tutorials/README.md): run your first simulation, choose a numerical method, or learn to develop with the ALPS libraries. It includes recommended learning paths, notebooks, and focused API examples.

## Installation

Build this checkout with standard CMake and Python packaging commands; see [CONTRIBUTING.md](CONTRIBUTING.md#getting-started-with-the-code) for dependencies and platform instructions. Use your existing dependency provider and toolchain. The [ALPS installation website](https://alps.comp-phys.org/install/) also covers upstream releases and Spack; its combined Boost.Python source-build instructions use a different interface from this SDK checkout.

Source builds require **CMake 3.27 or newer**. If your system provides an older version, follow the [CMake and Ninja setup](CONTRIBUTING.md#install-cmake-and-ninja) for pip installation or official binary downloads. For the full build instructions and exported CMake SDK target, see [the development setup](CONTRIBUTING.md#getting-started-with-the-code).

The [pyalps build instructions](python/pyalps/README.md) explain how to build the Python package against an installed C++ SDK.

See [CHANGELOG.md](CHANGELOG.md) for release notes and upgrade guidance.

## Develop from source

After installing the [prerequisites](CONTRIBUTING.md#prerequisites), configure, build and test the SDK:

```sh
cmake --preset default
cmake --build --preset default --parallel 2
ctest --preset default
cmake --install _build/default
```

Then follow the [editable Python installation](python/pyalps/README.md#editable-development), pointing `ALPS_DIR` at `_build/default/install/share/alps`.

The [module layout](src/alps/README.md) prepares shared libraries and applications for ALPSCore reconciliation. Utilities, HDF5, params, run configuration, Osiris, XML command-line parsing and modern ALEA statistics have separate exported libraries; both the MaxEnt solver and executable link components without `ALPS::alps`. `ALPS::headers` still provides a shared compile interface, so source ownership does not imply that every module is independent.

MaxEnt, segment CT-HYB and CT-INT use application-owned TOML schemas with separate scientific parameters, numerical input, output and execution settings. Their existing scientific methods remain; focused regression tests cover correctness fixes exposed by the migration. See the [Python configuration API](python/pyalps/README.md#typed-params-and-toml-migration). The remaining scheduler, parapack and DMFT workflows are still being migrated; their scientific model and lattice XML resources remain supported.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for development and contribution guidance.

## License and citations

ALPS is distributed under the terms in [LICENSE.txt](LICENSE.txt). See [CITATION.md](CITATION.md) for citation guidance.
