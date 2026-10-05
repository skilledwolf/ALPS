# Native ALEA analysis commands

With the pyalps bindings installed, the source-tree commands retain their names
and options:

```sh
python src/tools/alea/mean.py -v results.h5
python src/tools/alea/variance.py -v -n Energy results.h5
python src/tools/alea/variance.py -w -p /simulation/results results.h5
```

`-n/--name` selects logical observable names and may be repeated;
`-p/--path` selects the result group. Without `-n`, all results are selected.
`-v/--verbose` prints estimates and `-w/--write` stores them at `mean/value` or
`variance/value`. Inputs remain read-only unless `-w` is supplied. Unsupported
selected results are rejected before any estimates are written to that file.

Both commands use `pyalps.alea.read_result(archive, path)`, also used by
`loadMeasurements`. This preserves the concrete result type and its batches,
covariance, complex/elliptic uncertainties and observation weights. The commands
report that result's native `mean` or `variance`; they do not treat unequally
weighted bins as equally weighted samples. Variance describes the observations
represented by the estimator, which may themselves be batches. It does not
reconstruct individual observations within a bin. A mean-only result cannot
supply a variance.

Convert released files offline first, using the
[archive converter](../hdf5/README.md). `--alea-batches GROUP` recovers complete
ALPS linear histories as native batch results; `--core-alea KIND GROUP` converts
released ALPSCore ALEA results. Generic leaf conversion or `--alea` alone does
not create native statistical evidence from a legacy summary. Accumulator
checkpoints are not accepted as result groups.
