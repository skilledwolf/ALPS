# Copyright (C) 2015, 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
"""Run explicit native simplemc tasks; snapshots are emitted directly as VTK."""
import argparse
from pyalps.run_io import execute, write_run_files


def runs(group):
    temperatures = (5., 4.5, 4., 3.5, 3., 2.9, 2.8, 2.7, 2.6, 2.5, 2.4,
                    2.3, 2.2, 2.1, 2., 1.9, 1.8, 1.7, 1.6, 1.5, 1.2)
    cases = {
        "a": [("ising", "square lattice", length, temperature, 1., 0.)
              for length in (8, 16, 24) for temperature in temperatures],
        "b": [("ising", "square lattice", 128, temperature, 1., 0.)
              for temperature in (3., 2.3, 2.)],
        "c": [("xy", "square lattice", 64, .01, 1., 0.)],
        "d": [("heisenberg", "triangular lattice", 48, .01, -1., field)
              for field in (0., 5.)],
    }
    for index, (model, lattice, length, temperature, coupling, field) in enumerate(cases[group], 1):
        parameters = dict(ALGORITHM=model, LATTICE=lattice, L=length, T=temperature,
                          J=coupling, H=field, SWEEPS=65536 if group == "a" else 16384)
        output = {"results": f"parm9{group}.task{index}.out.h5"}
        execution = {}
        if group != "a":
            parameters["THERMALIZATION"] = 0
            output["snapshot_prefix"] = f"parm9{group}.task{index}"
            execution["snapshot_interval"] = 16384
        yield dict(parameters=parameters, output=output, execution=execution)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("group", choices=("a", "b", "c", "d"))
    arguments = parser.parse_args()
    job = write_run_files("parm9" + arguments.group, runs(arguments.group), baseseed=42)
    execute("simplemc", job)
