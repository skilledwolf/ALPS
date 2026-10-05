# Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
"""Follow running estimates by resuming complete native checkpoints.

Prefixes share samples, so their error bars are correlated. This is a visual
convergence diagnostic, not a stationarity test or a raw measurement series.
Repeat with longer warm-up and independent seeds to assess equilibration.
"""
import matplotlib.pyplot as plt
from pyalps import alea, hdf5
from pyalps.run_io import execute, write_run_file

parameters = {"LATTICE": "square lattice", "MODEL": "Ising", "UPDATE": "local",
              "L": 48, "J": [1.0], "T": 2.269186, "THERMALIZATION": 10000,
              "SWEEPS": 50000}
checkpoint = None
points = []
for segment in range(12):
    prefix = f"parm1a.segment{segment + 1}"
    run = write_run_file(prefix + ".toml", parameters=parameters,
        execution={"seed": 42, "max_sweeps": 5000, "bins": 128},
        input={} if checkpoint is None else {"checkpoint": checkpoint},
        output={"checkpoint": prefix + ".checkpoint.h5", "results": prefix + ".out.h5"},
        overwrite=True)
    filename = execute("spinmc", run)[0]
    checkpoint = prefix + ".checkpoint.h5"
    with hdf5.archive(filename) as archive:
        result = alea.BatchResult.read(archive, "/simulation/results/|Magnetization|")
    if result.count:
        points.append((result.count, result.mean[0], result.error[0]))

counts, means, errors = zip(*points)
plt.errorbar(counts, means, yerr=errors, marker="o")
plt.xlabel("Accumulated production updates")
plt.ylabel("Running mean of |Magnetization|")
plt.show()
