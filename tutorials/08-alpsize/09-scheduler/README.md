# Native scheduling and continuation: Wolff

`wolff` now uses the shared `alps::mc` runner and `alps::mcbase`, with the same
cluster update and Binder analysis as lessons 07 and 08. It accepts the
established TOML CLI:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/path/to/alps/install
cmake --build build --target wolff
build/wolff --validate run.toml
build/wolff run.toml
```

`--schema` prints the full schema. Use multiple run files for temperature scans.
The included run has 16 independent chains; `execution.chains` replaces
`NUM_CLONES`, and `execution.seed` is the base seed (chain IDs add their offset).
MT19937 and lagged Fibonacci RNGs are supported. An MPI-enabled SDK distributes
chains through the same native runner. Statistical pooling retains independent
partial batches; autocorrelation diagnostics remain separate per chain.

To interrupt reproducibly, set `execution.max_sweeps` to a per-chain budget,
including warmup. `execution.time_limit` provides a wall-clock limit. Output
checkpoints contain RNG state, every spin, neighbor ordering, sweep count, chain
identity and native measurement accumulators, including unfinished batches.
A warmup-only stop still writes a usable checkpoint and an analysis file with
an explicit unavailable-statistics reason.

To continue, create another run with `input.checkpoint` naming the checkpoint,
new output paths, and `execution.max_sweeps=0` (the default). Model, lattice,
seed, chain count, bins and warmup must agree. A completed run can be extended by
increasing `SWEEPS`; specify `THERMALIZATION` explicitly so its default does not
change with the target. Checkpoint validation completes before output publication.
Changed topology, invalid spins or inconsistent measurement counts are rejected.

Results preserve `<m^2>^2/<m^4>` as the Binder convention. `/simulation/joint`
retains the pooled batch evidence; per-chain diagnostics are under
`/simulation/realizations/0/clones/<id>/autocorrelation`. These native checkpoints
support exact continuation. Released Parapack/XDR checkpoints still require a
physical-state converter; they are not accepted as native checkpoints.

The separate `hello` example still demonstrates the old Parapack interface and
has not yet been migrated. Build the `wolff` target for this lesson.
