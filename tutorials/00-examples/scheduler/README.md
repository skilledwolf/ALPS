# Native Ising scheduling examples

Build the `ising1`, `ising2`, and `ising3` targets from the parent examples
project against an installed ALPS SDK. All three use consolidated parameters,
native ALEA statistics and the shared `alps::mc` runner:

```sh
build/scheduler/ising1 --schema
build/scheduler/ising1 --validate scheduler/run.toml
build/scheduler/ising1 scheduler/run.toml
build/scheduler/ising2 scheduler/lattice.toml
```

`ising1` and `ising3` simulate a periodic Ising chain. Their former scheduler
entry points now share one implementation and the same TOML CLI. `ising2` uses
the lattice library: specify exactly one `parameters.LATTICE` or
`parameters.GRAPH`, with optional `input.lattice_library`. Couplings are one,
spins are ±1, and the Hamiltonian is `-sum_b s_i*s_j`. Disordered lattices remain
unsupported. One sweep makes N random-site Metropolis proposals, with replacement.
Each proposal keeps or flips the selected spin with equal probability before
the Metropolis acceptance test. Allowing the current spin prevents an even
measurement stride from trapping the simulation in a parity sector when every
flip would be accepted (the infinite-temperature limit).

The original measurements are retained: energy and magnetization per site in
all variants; the chain variants also measure second and fourth magnetization
moments and `Correlations[d] = sum_i s_i*s_(i+d)/N`. A joint native batch result
under `/simulation/joint` retains covariance for later analysis; individual
results are under `/simulation/results`. Per-chain vector autocorrelation
hierarchies are under `/simulation/realizations/0/clones/<id>/autocorrelation/Moments`,
with components in the same order (energy, magnetization, then second/fourth
moments and all correlation distances for the chain variants).

Use multiple TOML files for scans. `execution.chains` selects independent
chains; the shared runner distributes them when built with MPI. `execution.rng`
selects `mt19937` or `lagged_fibonacci607`. Output paths are relative to the run
file. `execution.max_sweeps` sets a per-chain invocation budget including warmup;
`execution.time_limit` sets a wall-clock limit. Native checkpoints retain all
spins, ordered topology, RNG state, sweep counts and unfinished statistical
batches. Continue with `input.checkpoint`; increasing `SWEEPS` also extends a
completed run. Other physical parameters, warmup, chain count and statistical
configuration must agree. A warmup-only stop records unavailable statistics
and still writes a resumable checkpoint.

Old scheduler flags, parameter expressions and XML input are replaced by typed
TOML. Released scheduler/XDR physical checkpoints are not native checkpoints;
conversion and continuation support remains outstanding. The `evaluate` and `evaluate2` commands now build from one native postprocessor:

```sh
build/scheduler/evaluate --validate scheduler/evaluate.toml
build/scheduler/evaluate scheduler/evaluate.toml
```

`input.results` names a native chain results file from `ising1` or `ising3`;
`output.results` is a separate analysis file. Input is never modified. The
output retains joint batch evidence and includes correlations plus the original
`Binder cumulant of Magnetization` label. Despite that historical label, the
formulas are ratios: `evaluate` computes `<m^4>/<m^2>^2`; `evaluate2` computes
three times that ratio. They intentionally differ from the Wolff lesson.
Both use joint weighted jackknife propagation, including partial batches.
`--validate` checks the archive and propagation without writing. Insufficient
batches or nonpositive second moments in any jackknife sample reject analysis
before replacing output. Released separate-observable summaries do not establish
joint covariance and are not accepted as native joint evidence.
