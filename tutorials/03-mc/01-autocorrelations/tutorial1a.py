# Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
"""Compare native logarithmic blocking errors and autocorrelation estimates."""
import matplotlib.pyplot as plt
from pyalps import alea, hdf5
from pyalps.run_io import execute, write_run_files


def run(update="local", thermalization=10000, sweeps=50000, prefix="parm1a"):
    sizes = [2, 4, 8, 16, 32, 48]
    runs = [{"parameters": {"LATTICE": "square lattice", "L": size,
        "T": 2.269186, "J": [1.0], "THERMALIZATION": thermalization,
        "SWEEPS": sweeps, "UPDATE": update, "MODEL": "Ising"},
        "execution": {"seed": 42 + size},
        "output": {"results": f"{prefix}.L{size}.out.h5"}} for size in sizes]
    files = execute("spinmc", write_run_files(prefix, runs, overwrite=True))
    plt.figure()
    for size, filename in zip(sizes, files):
        with hdf5.archive(filename) as archive:
            result = alea.AutocorrelationResult.read(archive,
                "/simulation/realizations/0/clones/0/autocorrelation/|Magnetization|")
        # Coarse levels with few effective bins have unreliable uncertainty.
        levels = [i for i in range(result.levels) if result.level(i).observations >= 16]
        plt.plot([2**i for i in levels], [result.level(i).error[0] for i in levels],
                 marker="o", label=f"L={size}")
        print(f"L={size}: tau={result.tau if result.tau_available else 'unavailable'}")
    plt.xscale("log", base=2)
    plt.xlabel("Blocking width (production updates)")
    plt.ylabel("Native error of |Magnetization|")
    plt.legend()
    plt.show()


if __name__ == "__main__":
    run()
