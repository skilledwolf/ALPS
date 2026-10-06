# Quantum loop example

The former `loop_single` example program is replaced by the installed `loop`
application, which implements the same continuous-time loop algorithm and
observables. With the installation's `bin` directory on `PATH`, run:

```sh
loop --validate run.toml
loop run.toml
python runs.py scan
python runs.py disorder
```

The explicit spin model has spin 1/2 and antiferromagnetic isotropic exchange
`J = 1`. Energy and staggered magnetization squared are extensive. Energy
density and uniform/staggered susceptibilities use the geometric lattice
volume, including depleted sites. `Number of Sites` counts surviving spins.
The staggered susceptibility includes the imaginary-time integral.

The disorder example generates 100 separate TOML runs with recorded disorder
and depletion seeds. `execution.chains` runs independent Monte Carlo histories
of one quenched model; it does not generate additional disorder realizations.
Average each realization's result with equal weight for a quenched disorder
average, retaining the variation between realizations when estimating errors.

Native checkpoints resume through `input.checkpoint`. Resume released
Parapack checkpoints with ALPS 3.0 and convert their results for analysis.
