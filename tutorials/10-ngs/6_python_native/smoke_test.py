#!/usr/bin/env python3
"""Check Python Ising measurements and complete, transactional restart state."""

import copy
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import numpy as np

import ising
import pyalps.hdf5 as hdf5
from pyalps.alea import BatchResult
from pyalps.ngs import params


PARAMETERS = {"L": 7, "THERMALIZATION": 11, "SWEEPS": 701, "T": 2.0, "SEED": 42}
OBSERVABLES = ["Correlations", "Energy", "Magnetization",
               "Magnetization^2", "Magnetization^4"]
BASE = "/simulation/realizations/0/clones/0"


def assert_results_equal(actual, expected):
    assert sorted(actual) == sorted(expected) == OBSERVABLES
    for name in OBSERVABLES:
        assert actual[name].count == expected[name].count, name
        for field in ("mean", "error", "covariance", "batch_sums", "batch_counts"):
            np.testing.assert_array_equal(getattr(actual[name], field), getattr(expected[name], field))


simulation = ising.sim(PARAMETERS)
samples = {name: [] for name in OBSERVABLES}
assert simulation.fraction_completed() == 0.0
for _ in range(PARAMETERS["THERMALIZATION"] + PARAMETERS["SWEEPS"]):
    simulation.update()
    simulation.measure()
    if simulation.sweeps > PARAMETERS["THERMALIZATION"]:
        # Independently compute the physical observables from the spin stream.
        magnetization = simulation.spins.mean()
        values = {
            "Energy": -np.dot(simulation.spins, np.roll(simulation.spins, 1)) / simulation.length,
            "Magnetization": magnetization,
            "Magnetization^2": magnetization**2,
            "Magnetization^4": magnetization**4,
            "Correlations": [np.dot(simulation.spins, np.roll(simulation.spins, d)) / simulation.length
                             for d in range(simulation.length)],
        }
        for name, value in values.items():
            samples[name].append(np.atleast_1d(value))

results = simulation.collectResults()
assert simulation.fraction_completed() == 1.0
for name in OBSERVABLES:
    result = results[name]
    assert result.count == PARAMETERS["SWEEPS"], name
    components = PARAMETERS["L"] if name == "Correlations" else 1
    assert result.mean.shape == result.error.shape == (components,)
    np.testing.assert_allclose(result.mean, np.mean(samples[name], axis=0), rtol=0, atol=1e-14)
    # Independent weighted batch covariance and standard-error calculation.
    counts = result.batch_counts.astype(float)
    sums = result.batch_sums[counts > 0]
    counts = counts[counts > 0]
    mean = sums.sum(axis=0) / counts.sum()
    residuals = sums / counts[:, None] - mean
    count2 = np.dot(counts, counts)
    covariance = residuals.T @ (counts[:, None] * residuals) / (counts.sum() - count2 / counts.sum())
    error = np.sqrt(np.diag(covariance) * count2 / counts.sum()**2)
    np.testing.assert_allclose(result.covariance, covariance, rtol=1e-12, atol=1e-14)
    np.testing.assert_allclose(result.error, error, rtol=1e-12, atol=1e-14)

with tempfile.TemporaryDirectory() as directory:
    output = Path(directory) / "ising.h5"
    checkpoint = output.with_suffix(".clone0.h5")
    stopped = ising.sim(PARAMETERS)
    assert not stopped.run(lambda: stopped.sweeps == 148)
    assert stopped.fraction_completed() < 1.0
    hdf5.save_checkpoint(str(checkpoint), stopped.save)
    protocol = Path(directory) / "archive-protocol.h5"
    with hdf5.archive(protocol, "w") as archive:
        archive.set_context("/caller")
        archive["/"] = stopped
        assert archive.context == "/caller"

    with hdf5.archive(checkpoint, "r") as archive:
        assert archive["/parameters/format"] == "alps.params.v2"
        assert sorted(archive.list_children(BASE + "/measurements")) == OBSERVABLES
        assert sorted(archive.list_children(BASE + "/checkpoint")) == ["engine", "name", "spins", "sweeps"]
        # 137 samples overflow the 64 slots twice and leave an unfinished batch.
        assert archive[BASE + "/measurements/Energy/cursor/level"] >= 2
        counts = archive[BASE + "/measurements/Energy/batch/count"]
        target = 2**int(archive[BASE + "/measurements/Energy/cursor/level"])
        assert np.any((counts > 0) & (counts < target))
        assert archive[BASE + "/measurements/Energy/@kind"] == 6

    restored = ising.sim(dict(PARAMETERS, L=4, T=1.0, SEED=7))
    with hdf5.archive(protocol, "r") as archive:
        archive.set_context("/caller")
        restored.load(archive)
        assert archive.context == "/caller"
    assert restored.parameters == PARAMETERS
    assert restored.sweeps == stopped.sweeps
    assert restored.length == PARAMETERS["L"]
    np.testing.assert_array_equal(restored.spins, stopped.spins)
    assert_results_equal(restored.collectResults(), stopped.collectResults())
    assert restored.run(lambda: False)
    assert restored.sweeps == simulation.sweeps
    assert restored.fraction_completed() == simulation.fraction_completed()
    np.testing.assert_array_equal(restored.spins, simulation.spins)
    assert_results_equal(restored.collectResults(), results)
    original_rng, restored_rng = copy.deepcopy(simulation.random), copy.deepcopy(restored.random)
    assert [original_rng() for _ in range(8)] == [restored_rng() for _ in range(8)]
    assert restored.run(lambda: False)  # Already complete: no extra update/sample.
    assert_results_equal(restored.collectResults(), results)

    malformed = {
        BASE + "/checkpoint/spins": stopped.spins[:-1],
        BASE + "/checkpoint/sweeps": stopped.sweeps + 1,
        BASE + "/checkpoint/engine": "not a random engine",
        BASE + "/measurements/Energy/batch/sum": np.zeros((64, 2), dtype=np.float64),
    }
    for field, value in malformed.items():
        broken = Path(directory) / "broken.h5"
        shutil.copyfile(checkpoint, broken)
        with hdf5.archive(broken, "a") as archive:
            archive[field] = value
        previous = restored.__dict__
        previous_results = restored.collectResults()
        previous_rng = copy.deepcopy(restored.random)
        with hdf5.archive(broken, "r") as archive:
            archive.set_context("/caller")
            try:
                restored.load(archive)
            except (RuntimeError, ValueError, hdf5.ArchiveError):
                pass
            else:
                raise AssertionError(f"accepted malformed checkpoint field {field}")
            assert archive.context == "/caller"
        assert restored.__dict__ is previous
        assert restored.parameters == PARAMETERS
        assert restored.sweeps == simulation.sweeps
        np.testing.assert_array_equal(restored.spins, simulation.spins)
        assert_results_equal(restored.collectResults(), previous_results)
        current_rng = copy.deepcopy(restored.random)
        assert [current_rng() for _ in range(8)] == [previous_rng() for _ in range(8)]

    # Exercise the actual -c CLI and the result consumer using restored params.
    subprocess.run([sys.executable, str(Path(__file__).with_name("main.py")), "-c", str(output)],
                   check=True, capture_output=True, text=True)
    with hdf5.archive(output, "r") as archive:
        with archive.native() as native:
            assert params(native) == PARAMETERS
            loaded_results = {name: BatchResult.read(native, "/simulation/results/" + name)
                              for name in OBSERVABLES}
        assert_results_equal(loaded_results, results)
        for name, result in results.items():
            path = "/simulation/results/" + name
            assert archive[path + "/@kind"] == 5
            np.testing.assert_array_equal(archive[path + "/mean/value"], result.mean)
            np.testing.assert_array_equal(archive[path + "/mean/error"], result.error)
            np.testing.assert_array_equal(archive[path + "/batch/sum"], result.batch_sums)
            np.testing.assert_array_equal(archive[path + "/batch/count"], result.batch_counts)

print("native python simulation: ok")
