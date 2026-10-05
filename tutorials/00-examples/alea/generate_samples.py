"""Generate chronological AR(1) observations for the Python analysis examples."""
import numpy as np
from pyalps import hdf5

rng = np.random.default_rng(42)
samples = np.zeros((32768 + 1024, 2))
for t in range(1, len(samples)):
    samples[t] = 0.95 * samples[t-1] + rng.normal(size=2)
# Discard startup relaxation; preserve every subsequent observation in order.
with hdf5.archive("timeseries.h5", "w") as archive:
    archive["/samples/E"] = samples[1024:, 0]
    archive["/samples/m"] = samples[1024:, 1]
