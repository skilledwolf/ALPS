# Library examples

These small programs demonstrate individual ALPS APIs. For complete simulation workflows and a recommended learning sequence, start with the [tutorial guide](../README.md).

| Directory | Subject |
| --- | --- |
| [alea](alea/) | Statistical estimates and error analysis |
| [hdf5](hdf5/) | Saving custom C++ types in HDF5 archives |
| [ietl](ietl/) | Iterative eigenvalue methods |
| [model](model/) | Symbolic and numerical model construction |
| [parapack](parapack/) | Simulation workers and parameter scans |
| [sampling](sampling/) | Sampling distributions |
| [scheduler](scheduler/) | Scheduled simulations |
| [fortran](fortran/) | Calling Fortran simulation code through the C++ bridge |

## Build against an installed SDK

These examples require CMake 3.27 or newer; see the [CMake setup instructions](../../CONTRIBUTING.md#install-cmake-and-ninja) if your system provides an older version.

From this directory, build the C++ examples together:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/path/to/alps
cmake --build build --parallel 2
ctest --test-dir build --output-on-failure
```

Fortran examples are a separate project and require a Fortran compiler:

```sh
cmake -S fortran -B build-fortran -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/path/to/alps
cmake --build build-fortran --parallel 2
ctest --test-dir build-fortran --output-on-failure
```

Both projects register tests by default; pass `-DALPS_BUILD_TESTING=OFF` to build without registering tests. These standalone builds replace the root `ALPS_BUILD_EXAMPLES` option. To use a local ALPS checkout, build and install its SDK first, then pass that installation prefix in the commands above.

CTest locates the SDK's XML resources. When running an executable directly, use the input files in its example directory; set `ALPS_ROOT` to the SDK prefix if the SDK has been relocated.

## Native classical examples

`parapack/ising/ising` and `parapack/heisenberg/heisenberg` keep their command
names and use the shared TOML CLI: `program --schema`, `program --validate
run.toml`, or `program run.toml`. Each directory contains a runnable replica
exchange example. Both now use native ALEA and the shared replica-exchange
controller; their old Parapack worker/evaluator classes are removed.

Choose a lattice or graph, a default bond coupling `J`, and optional bond-type
couplings `J0`, `J1`, etc. The energy is `-sum_b J_b s_i*s_j` (the dot product
for Heisenberg unit vectors). Ising uses heat-bath updates and Heisenberg uses
isotropic Metropolis proposals. Self-bonds contribute constant energy. Energy
is extensive; magnetization moments are per site. Heisenberg retains both the
vector-magnitude and z-component second/fourth moments and Binder ratios.
Centered energy moments supply heat capacity and all nonlinear estimates use
joint batches.

For fixed temperature set `T` or `BETA`. Omitting both selects the classical
infinite-temperature limit, `BETA = 0`. For exchange set `ALGORITHM` to
`"ising; exchange"` or `"heisenberg; exchange"`, then supply explicit
`TEMPERATURE_SET`/`INVERSE_TEMPERATURE_SET` arrays or `NUM_REPLICAS` and
`T_MIN`/`T_MAX` or `BETA_MIN`/`BETA_MAX`. The shared schema lists exchange
intervals, randomized ordering, disabled exchange, and rate/population feedback.
Adaptive ladders require `execution.chains = 1`. Uniform-temperature grids and
population feedback require positive inverse temperatures; explicit and
inverse-temperature grids also allow the classical zero endpoint. At zero beta,
results report `BETA = 0` and omit infinite-valued temperature diagnostics.

Native checkpoints retain all physical walkers, RNG streams, partial batches
and feedback history. Production can be extended with a larger `SWEEPS`.
MPI distributes independent chains/ladders and permits changing process count
on restart. Set `execution.parallel = "replicas"` to advance physical walkers concurrently
across MPI ranks; `"chains"` remains the default independent-ladder layout.
Each temperature keeps its original sample order and autocorrelation history,
and checkpoints can resume in either layout at a different rank count.
Ranks store only their owned physical walkers; native checkpoint groups stream
through a bounded buffer into root's disk cache. Temperature histories remain
replicated to preserve their chronology. The separate legacy `exchange` example
still supplies nested replica/spatial execution until that port is complete.

The [quantum loop example](parapack/loop/README.md) runs the installed `loop`
application, which replaces the former `loop_single` example program. Its
explicit disorder jobs retain separate quenched realizations.

The [classical energy Wang–Landau example](parapack/wanglandau/README.md) also
uses native TOML runs and ALEA. It retains density-of-states learning, overlapping
window stitching, fixed-weight microcanonical sampling, temperature reweighting,
reference-normalized entropy/free energy, histograms and exact continuation.
