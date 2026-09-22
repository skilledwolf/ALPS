# Develop ALPS with binary dependencies

Clone the repository, then use one build command. The managed setup downloads dependencies, builds and installs the C++ SDK, and installs editable Python bindings. Only ALPS and its bindings are compiled. Existing manual presets remain available for MPI/HPC installations and custom toolchains.

## Linux and macOS

Install [Pixi](https://pixi.sh/latest/installation/) once. On macOS, install Apple's Command Line Tools once with `xcode-select --install` to provide the platform SDK. Then, from the checkout:

```sh
pixi run --locked dev
```

The committed `pixi.lock` supplies binary Python, CMake, Ninja, compilers, Boost, serial HDF5, OpenBLAS/LAPACK, and Python dependencies. Linux x64/ARM64 and macOS Intel/ARM64 are included. The lockfile targets Linux with glibc 2.28+/kernel 4.18+ and macOS 13+. Linux uses an environment-provided compiler and sysroot; system development packages are unnecessary. The default build disables MPI.

## Windows

Install standard, GIL-enabled 64-bit CPython 3.12+ and Visual Studio's **Desktop development with C++** workload once. Include the target tools for your architecture. From an ordinary PowerShell in the checkout:

```powershell
python tools/dev.py
```

The helper finds and activates MSVC, creates a checkout-local virtual environment, and downloads the checksum-pinned dependency SDK used by CI. It does not clone or bootstrap vcpkg, compile dependencies, or change your global Python environment. The selected architecture follows Python: x64 Python builds x64 ALPS, including on an ARM64 computer. Use ARM64 Python for a native ARM64 build.

Windows reuses the existing binary exports because the needed conda-forge packages do not cover Windows ARM64. The Unix and Windows paths share the same build/development helper.

## Everyday commands

| Operation | Linux/macOS | Windows |
| --- | --- | --- |
| Build SDK and editable bindings | `pixi run --locked dev` | `python tools/dev.py` |
| Build and run native/Python/SDK tests | `pixi run --locked test` | `python tools/dev.py test` |
| Configure without compiling ALPS | `pixi run --locked python tools/dev.py configure` | `python tools/dev.py configure` |
| Build without Python bindings | `pixi run --locked python tools/dev.py --cpp-only` | `python tools/dev.py --cpp-only` |
| Run with the installed SDK programs on PATH | `pixi run --locked python tools/dev.py run <command>` | `python tools/dev.py run <command>` |

For example, replace `<command>` with `python my_simulation.py`, `spinmc --help`, or `python -c "import pyalps; print(pyalps.__file__)"`. The `run` command uses the existing build and does not rebuild it. On Unix, `pixi shell` also exposes the Python environment for interactive work.

Python source edits take effect in a new interpreter. After C++ changes, rerun the development command: it rebuilds incrementally, reinstalls the SDK, and refreshes the bindings in a persistent build directory. The ordinary setup builds the SDK, applications, and bindings; `test` additionally builds the regression executables. The default is two compiler processes and two concurrent tests, suitable for smaller laptops. Put `--jobs 4` before `test` or `configure` to increase compilation parallelism on larger machines. After intentionally changing compilers or toolchains, use `--fresh` before the command to reset the SDK and binding CMake configurations. Compiler caching is enabled if `ccache` is on PATH; Pixi includes it.

Configure your editor to use `.pixi/envs/default/bin/python` on Unix or `_build/dev/windows-<architecture>/venv/Scripts/python.exe` on Windows. Compilation databases are in `_build/dev/<platform>/sdk/compile_commands.json`; the installed SDK is under the neighboring `install/` directory. The helper prints that path. Generated files stay under `_build/` and `.pixi/`, independently of existing manual presets and local virtual environments.

## Maintain the environments

Use `pixi update` to intentionally refresh the Unix lockfile, then validate the developer setup on Linux and macOS. Ordinary setup uses `--locked` so it cannot silently change dependencies. Windows Python requirements are in `tools/dev-requirements.txt` and are installed as binary wheels. The helper updates its isolated environment when that file changes. Native Windows library pins come from `.github/dependencies.json`; changing a checksum selects a new verified dependency directory. Changes to the lockfile or Windows dependency inputs automatically refresh the SDK and binding CMake configurations, preventing stale package locations.

The binary exports currently come from the repository recorded in that manifest. Before moving upstream ownership, publish a complete dependency release under `ALPSim/ALPS`, import its checksums, and validate it. Do not change the repository field to a location that has not published the archives. See [CI dependency maintenance](../CONTRIBUTING.md#ci-dependency-binaries).
