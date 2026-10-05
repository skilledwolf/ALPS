# Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
"""Compare native blocked uncertainties for the same Monte Carlo trajectory.

Each requested batch capacity repeats the same seeded trajectory. Result bins
are weighted, adaptively merged blocks; their slots do not encode time order.
The x axis reports the actual average occupied-block width, including partial
blocks, rather than claiming access to a historical binning hierarchy.
"""
import matplotlib.pyplot as plt
from pyalps import alea, hdf5
from pyalps.run_io import execute, write_run_files


def run(update="local", thermalization=10000, sweeps=50000, prefix="parm1a"):
    runs = []
    for size in [2, 4, 8, 16, 32, 48]:
        for bins in [64, 128, 256, 512, 1024]:
            runs.append({"parameters": {"LATTICE": "square lattice", "L": size,
                "T": 2.269186, "J": [1.0], "THERMALIZATION": thermalization,
                "SWEEPS": sweeps, "UPDATE": update, "MODEL": "Ising"},
                # Accumulator capacity does not affect updates or RNG draws.
                "execution": {"seed": 42 + size, "bins": bins},
                "output": {"results": f"{prefix}.L{size}.bins{bins}.out.h5"}})
    manifest = write_run_files(prefix, runs, overwrite=True)
    files = execute("spinmc", manifest)
    plt.figure()
    for size in [2, 4, 8, 16, 32, 48]:
        points = []
        for run_config, filename in zip(runs, files):
            if run_config["parameters"]["L"] != size:
                continue
            with hdf5.archive(filename) as archive:
                result = alea.BatchResult.read(archive, "/simulation/results/|Magnetization|")
            occupied = result.batch_counts[result.batch_counts > 0]
            points.append((occupied.mean(), result.error[0]))
        points.sort()
        plt.plot(*zip(*points), marker="o", label=f"L={size}")
    plt.xscale("log", base=2)
    plt.xlabel("Mean occupied batch width (production updates)")
    plt.ylabel("Native error of |Magnetization|")
    plt.legend()
    plt.show()


if __name__ == "__main__":
    run()
