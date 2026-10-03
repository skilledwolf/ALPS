# Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
import pathlib
import subprocess
import sys
import tempfile

executable = sys.argv[1]
with tempfile.TemporaryDirectory(prefix="alps-cthyb-cli-") as temporary:
    root = pathlib.Path(temporary)
    (root / "delta.dat").write_text("0 -0.5 -0.5\n1 -0.5 -0.5\n2 -0.5 -0.5\n")
    document = """[parameters]
BETA = 2.0
N_ORBITALS = 2
N_TAU = 2
N_MEAS = 1
THERMALIZATION = 0
SWEEPS = 16
U = 1.0
MEASURE_time = false
[input]
delta = "delta.dat"
[output]
results = "result.h5"
"""
    run = root / "run.toml"
    run.write_text(document)
    work = root / "work"
    work.mkdir()

    def invoke(*arguments):
        return subprocess.run([executable, *map(str, arguments)], cwd=work,
                              capture_output=True, text=True, timeout=45)

    result = invoke("--validate", run)
    assert result.returncode == 0, result.stderr
    assert not (root / "result.h5").exists(), "Validation must not start a simulation"
    result = invoke(run)
    assert result.returncode == 0, result.stderr
    assert (root / "result.h5").stat().st_size > 0
    assert (root / "delta.dat").exists(), "Input file disappeared"
    assert not list(work.iterdir()), "Text output must be opt-in"
    run.write_text(document.replace("SWEEPS = 16", "SWEEPS = 0"))
    result = invoke("--validate", run)
    assert result.returncode != 0 and "SWEEPS" in result.stderr
    result = invoke("--time-limit", "5", run)
    assert result.returncode != 0 and "Unknown option" in result.stderr
    result = invoke(root / "result.h5")
    assert result.returncode != 0, "Legacy HDF5 parameter input must be rejected"
