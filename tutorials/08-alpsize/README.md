# Native Wolff lessons

Lessons `07-alea` and `08-lattice` demonstrate the same zero-field ferromagnetic
Ising Wolff update, with Hamiltonian `-sum_b s_i*s_j` and unit bond couplings.
The first constructs a periodic square lattice explicitly; the second uses the
ALPS lattice library. `wolff.hpp` shares their update, native ALEA analysis and
output code, so the lattice lesson does not duplicate the simulation.

Build either lesson with an installed SDK:

```sh
cmake -S 07-alea -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/path/to/alps/install
cmake --build build
build/wolff --validate 07-alea/run.toml
build/wolff 07-alea/run.toml
```

Use `08-lattice` instead of `07-alea` for the lattice-library simulation. Each
executable remains named `wolff`; inputs are now TOML files instead of parameter
text on standard input. `--schema` prints its schema and `--validate` resolves
parameters and constructs the lattice without sampling or writing output.
Output paths are relative to the run file. The default run has `L=32`, `T=2.2`,
32768 production cluster updates, seed 93812, and 64 statistical batches.
Warmup defaults to `SWEEPS/8` unless `THERMALIZATION` is supplied. A sweep here
means one cluster flip, not one attempted update per lattice site.

The lattice lesson accepts exactly one `parameters.LATTICE` or `parameters.GRAPH`,
additional lattice parameters, and an optional `input.lattice_library` path.
It constructs the lattice directly from the typed parameters. The run
configuration and saved parameters use the consolidated typed API. The separate
`lattice` executable in lesson 08 remains the original graph-iteration example.

Results are written atomically to `output.results`. The three magnetization
moments and `Binder Ratio of Magnetization` are native results under
`/simulation/results`. The Binder convention remains `<m^2>^2/<m^4>`; its error
uses the joint batch evidence stored at `/simulation/joint`, including any
partial batch. If a jackknife sample has zero fourth moment, the file records an
explicit reason under `/simulation/unavailable` rather than a spurious ratio.
Per-observable native autocorrelation levels are retained for
`pyalps.loadBinningAnalysis`; measurements load with `pyalps.loadMeasurements`.
Printed errors use batches; tau, when available, comes from the autocorrelation
hierarchy. These output files are analysis results, not physical restart states.

The obsolete `.ip`/`.op` pairs have been replaced by `run.toml` and physics/error
regression checks. Released statistical archives can be migrated with the
appropriate offline `alps-hdf5-convert` ALEA profile.

Lesson [09-scheduler](09-scheduler/README.md) adds native independent-chain
scheduling and exact physical/statistical checkpoints to the same Wolff update.
Its `hello` command introduces the typed TOML CLI without simulation machinery.
