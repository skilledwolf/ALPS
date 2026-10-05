# Python Ising simulation with native ALEA

The lesson implements random-site Metropolis sweeps of a periodic square Ising
lattice, with Hamiltonian `-sum_<ij> s_i s_j`. Energy and magnetization are per
site. Implement `step` and `observables` in `ising-skeleton.py`; the parameter,
RNG, statistics and output code is shared with `solution/ising.py`.

From this directory, with pyalps installed:

```sh
python solution/ising.py --schema
python solution/ising.py --validate run.toml
python solution/ising.py run.toml
python solution/ising_binder.py run.toml
```

The Binder script plots the saved ratio. `solution/run.py` remains an alternative
entry point. With no arguments, the solution scripts retain their original beta
scans; the Binder script scans sizes 4, 6 and 8. TOML files use the shared typed
parameter validator. Each file describes one run, with output relative to that
file and seed/bin settings in `[execution]`. Multiple files must have distinct
outputs and cannot overwrite their run files. Every instance owns its RNG.

One five-component `BatchAccumulator` records `[E, m, |m|, m^2, m^4]` in aligned
batches. Per-observable `AutocorrelationAccumulator`s retain the logarithmic binning
hierarchies. Both exclude thermalization. The printed autocorrelation time is
marked unavailable until the estimator has sufficient bins.

The Binder **ratio** is `<m^4>/<m^2>^2`, not the cumulant
`1 - <m^4>/(3 <m^2>^2)`. Its bias correction and uncertainty use joint jackknife
samples, retaining covariance between the two moments. It is computed before
plotting; dividing independent `FloatWithError` summaries would lose that
covariance. If any required jackknife denominator is zero, the result is omitted
and the reason is saved under `/simulation/unavailable`.

Files contain typed parameters and run settings, individual native results under
`/simulation/results`, the original joint batches under `/simulation/joint`, and
the autocorrelation hierarchies under
`/simulation/realizations/0/clones/0/autocorrelation`.
`pyalps.loadMeasurements` reads the named results and `loadBinningAnalysis` reads
the per-observable diagnostic curves; `pyalps.alea.read_result` reads
the joint result and diagnostics without losing their statistical evidence.
Saving uses atomic publication. These files are analysis output, not physical
restart checkpoints; this introductory lesson has no resume operation.
