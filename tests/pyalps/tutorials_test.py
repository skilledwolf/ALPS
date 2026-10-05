# Copyright (C) 2026 by the ALPS collaboration
# Part of the ALPS Project — see LICENSE.txt for full license text.
# SPDX-License-Identifier: MIT

"""Run the pure-Python ngs tutorials as part of the wheel test suite.

tutorials/10-ngs/5_export_python is covered by build.yml, which compiles it as a
downstream CMake project. Its two pure-Python siblings needed no compiler and so
had no coverage at all -- and both were broken in ways only running them shows:
6_python_native could not store an ngs.params or its measurements dict through
`archive[path] = ...`, and neither one's load() had ever executed (one wrote to
the archive instead of reading from it, both called double()).

Running them here rather than in build.yml is deliberate: they need an installed
pyalps, which is what this suite already has, so cibuildwheel exercises them
against every wheel it builds.
"""

import os
import runpy
import subprocess
import sys
from itertools import product
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest


TUTORIALS = Path(__file__).resolve().parents[2] / "tutorials" / "10-ngs"

# 5_export_python is absent: its smoke test imports a compiled extension that
# only exists after the downstream CMake build in build.yml.
PURE_PYTHON_TUTORIALS = ["6_python_native", "7_python_extend"]


@pytest.mark.parametrize("tutorial", [
    "10-ngs/6_python_native", "10-ngs/7_python_extend",
])
@pytest.mark.parametrize("length", [3, 4])
def test_ising_metropolis_uses_measured_energy(tutorial, length, monkeypatch):
    # Fresh globals isolate the distinct ising modules, including mcbase subclasses.
    update = runpy.run_path(str(TUTORIALS.parent / tutorial / "ising.py"))["sim"].update
    # Exercise one proposal of the actual kernel without constructing or reseeding a simulation.
    monkeypatch.setitem(update.__globals__, "range", lambda _: (0,))
    beta = 0.5
    for configuration in product((-1, 1), repeat=length):
        spins = np.array(configuration)
        energy = -np.dot(spins, np.roll(spins, -1))
        for site in range(length):
            flipped = spins.copy()
            flipped[site] *= -1
            delta = -np.dot(flipped, np.roll(flipped, -1)) - energy
            probability = min(1., np.exp(-beta * delta))
            uniforms = (np.nextafter(probability, 0.),
                        np.nextafter(probability, 1.)) if probability < 1. else (0., 0.75)
            for uniform in uniforms:
                random = iter(((site + 0.5) / length, uniform)).__next__
                state = SimpleNamespace(length=length, beta=beta, spins=spins.copy(), random=random)
                update(state)
                np.testing.assert_array_equal(
                    state.spins, flipped if uniform < probability else spins,
                    err_msg=f"{tutorial}: spins={configuration}, site={site}, delta={delta}, u={uniform}",
                )


@pytest.mark.parametrize("tutorial", PURE_PYTHON_TUTORIALS)
def test_tutorial_smoke_test(tutorial, tmp_path):
    directory = TUTORIALS / tutorial
    script = directory / "smoke_test.py"
    if not script.is_file():
        pytest.skip(f"{script} is not present in this tree")

    # Run from a scratch directory so nothing lands in the source tree, with
    # the tutorial itself importable (each has its own ising.py).
    environment = dict(os.environ)
    environment["PYTHONPATH"] = os.pathsep.join(
        [str(directory)] + ([environment["PYTHONPATH"]]
                            if environment.get("PYTHONPATH") else []))

    completed = subprocess.run(
        [sys.executable, str(script)],
        cwd=str(tmp_path),
        env=environment,
        capture_output=True,
        text=True,
        timeout=900,
    )

    assert completed.returncode == 0, (
        f"{tutorial}/smoke_test.py failed\n"
        f"--- stdout ---\n{completed.stdout}\n"
        f"--- stderr ---\n{completed.stderr}")
    assert completed.stdout.strip().endswith("ok"), completed.stdout


def test_every_pure_python_tutorial_has_a_smoke_test():
    """A new pure-Python tutorial should not be able to arrive untested."""
    if not TUTORIALS.is_dir():
        pytest.skip("tutorials/10-ngs is not present in this tree")

    untested = []
    for directory in sorted(TUTORIALS.iterdir()):
        if not directory.is_dir() or not list(directory.glob("*.py")):
            continue
        if list(directory.glob("*.cpp")):
            continue        # compiled tutorials are covered by build.yml
        if not (directory / "smoke_test.py").is_file():
            untested.append(directory.name)

    assert not untested, f"pure-Python tutorials without a smoke test: {untested}"
