# Classical energy Wang–Landau

`wanglandau` keeps its command name and uses the native TOML CLI:

```sh
wanglandau --schema
wanglandau --validate learn.toml
wanglandau learn.toml
wanglandau measure.toml
wanglandau reweight.toml
```

Run the three phases separately: validation reads each phase's input before
execution. The supplied example samples the nonpositive-energy window of an
8×8 square Ising lattice. Its reweighted results are restricted to that window;
include the full spectrum when the omitted positive energies matter.

The Hamiltonian is `E = -COUPLING * sum_b s_i*s_j`, with integer `COUPLING`
(including negative and zero), Ising spins ±1 and the selected lattice-library
graph. Self-bonds contribute constant energy and never a flip field. A symmetric
keep/flip proposal avoids an even-sweep parity trap at flat weights.

## Learning and windows

Set `MODE = "learn"` and an inclusive two-integer `ENERGY_RANGE`. Separate
`ENERGY_WALK_RANGE` and `ENERGY_MEASURE_RANGE` override it; the measurement range
must lie inside the walk range. Every attempted site update records its visited
energy. Outside the measurement window, the log-weight increment is multiplied
by `VISIT_PENALTY` (default 8). Only weights inside the measurement window are
exported for subsequent use.

`CHECK_INTERVAL`, `INITIAL_UPDATE_FACTOR`, `FINAL_UPDATE_FACTOR` and
`FLATNESS_THRESHOLD` control refinement. A flat histogram halves `log(f)`;
learning finishes after a flat stage at or below the final factor. Both endpoints
participate. Energy spectra can contain holes: flatness considers all energies
visited so far in the measurement range, and includes zero counts when a
previously discovered energy is absent in a later stage. As with any sampled
DOS, completion cannot prove that a rare energy sector has been discovered;
check independent runs and increase the checking interval when needed.

`input.weights = ["window1.h5", "window2.h5"]` in measurement mode stitches
completed, overlapping learning windows. Paths resolve relative to the TOML
file. Windows must describe the same Hamiltonian and agree on accessible
energies in overlaps; at least one visited energy anchors every overlap. The
mean log-weight difference aligns arbitrary additive offsets, including zero
and negative values. Disconnected or insufficient coverage rejects before
outputs are touched. An explicit visited mask distinguishes spectrum holes from
valid zero weights.

## Measurements and reweighting

`MODE = "measure"` freezes the input weights. `THERMALIZATION` sweeps precede
`SWEEPS` production samples. Each sample retains the joint vector of energy-bin
occupancy and normalized magnetization's first, second and fourth moments.
Zeros are retained for bins not visited on that sweep; joint batches therefore
preserve both cross-bin correlations and unequal bin occupancy. Walk-only bins
may use flat weights and contribute no microcanonical measurement.

`MODE = "reweight"` reads `input.measurements = ["microcanonical.h5"]` and a
positive `TEMPERATURE_SET`. Several measurement files may be pooled if they
used the same weights (up to an additive constant), Hamiltonian and measurement
range. Only complete measurements are accepted. Reweighting is a deterministic
analysis with one execution chain and no simulation checkpoint.

The importance weight is `exp(log(g(E)) - beta*E)`, multiplied by the measured
bin occupancy. Native joint jackknife propagation provides errors for energy,
heat capacity, magnetization moments and Binder ratio. Stable exponential scaling
and energy differences avoid large absolute offsets in partition sums. To
obtain absolute free energy and entropy, supply both `REFERENCE_BIN` and its
known `REFERENCE_LOGG = log(degeneracy)`. Without this normalization these four
absolute quantities are marked unavailable; normalized observables still work.
Bins absent from a resample also produce explicit unavailable estimates rather
than invented uncertainties.

## Native evidence and continuation

`/weights` contains the energy range, Hamiltonian identity, log weights, visited
mask and completion flag. Measurement files keep the native ALEA joint batch
result under `/microcanonical/joint`. Per-chain final/overall/ascending and
measurement histograms remain under `/simulation/realizations/0/clones/<id>`;
round-trip diagnostics are ordinary native results. Microcanonical and canonical
results use the existing `/simulation/replicas/<id>` layout and load with
`pyalps.loadMeasurements`.

`output.checkpoint` and `input.checkpoint` retain spins, RNG, refinement stage,
weights, histograms and unfinished batches. `execution.max_sweeps` or the time
limit can interrupt either sampling phase. Increasing `SWEEPS` extends fixed-
weight production. MPI distributes independent learning/measurement chains;
changing MPI rank count on restart preserves the evidence. It does not split a
single window's updates across ranks.

The old Parapack Wang–Landau worker, XDR intermediates and legacy observable
adapter have been removed. Finish released simulations and their intermediate
learning/measurement workflows with ALPS 3.0. Released final HDF5 observables
remain eligible for results-only offline conversion; converted analysis results
are not native simulation checkpoints.
