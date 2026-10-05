# Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
"""Inspect each chain's ordered bins and weighted running mean.

Bin averages smooth within-bin fluctuations. This is a visual diagnostic;
repeat with longer warm-up and independent seeds to assess equilibration.
"""
import matplotlib.pyplot as plt
import numpy as np
from pyalps import alea, hdf5
from pyalps.run_io import execute, write_run_file

parameters = {"LATTICE": "square lattice", "MODEL": "Ising", "UPDATE": "local",
              "L": 48, "J": [1.0], "T": 2.269186, "THERMALIZATION": 10000,
              "SWEEPS": 50000}
run = write_run_file("parm1a.toml", parameters=parameters,
    execution={"seed": 42, "bins": 128}, output={"results": "parm1a.out.h5"},
    overwrite=True)
filename = execute("spinmc", run)[0]
with hdf5.archive(filename) as archive:
    series = alea.BatchAccumulator.read(archive,
        "/simulation/realizations/0/clones/0/series/|Magnetization|")
result = series.result()
order = np.argsort(series.batch_offsets)
order = order[result.batch_counts[order] > 0]
counts = result.batch_counts[order]
sums = result.batch_sums[order, 0]
ends = series.batch_offsets[order] + counts
plt.plot(ends - counts / 2, sums / counts, label="Bin averages")
plt.plot(ends, np.cumsum(sums) / np.cumsum(counts), label="Running mean")
plt.xlabel("Production updates")
plt.ylabel("|Magnetization|")
plt.legend()
plt.show()
