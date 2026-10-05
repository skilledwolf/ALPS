# ALEA analysis examples

Run `python generate_samples.py` with PyALPS installed to create `timeseries.h5`
in the working directory. It contains raw chronological observations of two
stationary AR(1) processes, labelled `E` and `m` for these lessons. No physical
Ising interpretation is intended. Compressed bin means cannot substitute for
raw observations in lag-correlation or running-mean calculations.

Build the C++ examples through the parent directory with an installed SDK:

```sh
cmake -S .. -B build -DCMAKE_PREFIX_PATH=/path/to/alps/install
cmake --build build --target example_mean example_variance example_autocorrelation example_error example_running_mean
```

Run `build/alea/example_mean` (and the other four executables) in the directory
containing `timeseries.h5`. The corresponding Python programs use the same data.

- Mean and variance use native ALEA moments; variance has the N-1 correction.
- Running means use forward and reverse sums of ordinary arrays.
- Lag correlations use a zero-padded FFT, normalized by `(N-lag)*sample_variance`.
  Lag zero is excluded. The exponential fit uses positive correlations between
  the first 80% and 20% crossings of the lag-one value, excluding the lower
  crossing. Its amplitude refers to lag zero.
- The integrated time sums through the 20% crossing and adds the continuum
  exponential tail from the last included lag plus one half. The error example
  compares independent-sample, native ALEA binning, and fitted-correlation errors.

The shared `analysis.hpp` contains the small array calculations used in the
lessons; it introduces no statistical container or archive compatibility layer.
Mean, variance and corrected error are written below `/analysis` in the input
file. Raw `/samples` remain unchanged. Binning estimates depend on having enough
samples per level; a returned error does not establish convergence.
