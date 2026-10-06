# Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
"""Native spinmc must retain complete chains and aligned analysis evidence."""
import subprocess

import h5py
import numpy as np
import pytest

import pyalps
from pyalps import alea, hdf5
from pyalps.run_io import execute, write_run_file
from conftest import alps_program


@pytest.fixture
def executable():
    return alps_program("spinmc")


def run_file(directory, name, *, parameters=None, execution=None, input=None, output=None):
    p = dict(MODEL="Ising", LATTICE="chain lattice", L=5, T=1.8,
             J=[.8], SWEEPS=37, THERMALIZATION=3, UPDATE="local")
    p.update(parameters or {})
    if "GRAPH" in p:
        p.pop("LATTICE")
    e = dict(seed=2873, chains=2, bins=8)
    e.update(execution or {})
    o = dict(results=name + ".h5")
    o.update(output or {})
    return write_run_file(directory / (name + ".toml"), parameters=p,
                          execution=e, input=input, output=o)


def results(filename):
    with hdf5.archive(filename) as archive:
        return {pyalps.hdf5_name_decode(name): alea.BatchResult.read(
                    archive, "/simulation/results/" + name)
                for name in archive.list_children("/simulation/results")}


@pytest.mark.parametrize("parameters", [
    {"h": [-.17]},
    {"MODEL": "XY", "J": [.8, .6, -.2, .9], "D": [.2, .1], "h": [.1, -.2]},
    {"MODEL": "Heisenberg", "J": [.8, .1, .2, -.3, .6, .1, 0., .2, .9],
     "D": [.1, .2, .3], "h": [.1, 0., -.2]},
    {"MODEL": "O(4)", "UPDATE": "cluster"},
    {"UPDATE": "cluster", "J": [-.8], "L": 4},
    *[{"MODEL": "Potts", "q": q, "UPDATE": "cluster"} for q in (3, 4, 10)],
])
@pytest.mark.parametrize("rng", ["mt19937", "lagged_fibonacci607"])
def test_spinmc_restart_retains_every_native_result(executable, tmp_path, parameters, rng):
    full = run_file(tmp_path, "full", parameters=parameters,
                    execution={"rng": rng}, output={"checkpoint": "full-checkpoint.h5"})
    stopped = run_file(tmp_path, "stopped", parameters=parameters,
                       execution={"max_sweeps": 14, "rng": rng}, output={"checkpoint": "partial.h5"})
    execute(executable, [full, stopped])
    resumed = run_file(tmp_path, "resumed", parameters=parameters,
                       execution={"rng": rng}, input={"checkpoint": "partial.h5"})
    execute(executable, resumed)
    expected, actual = results(tmp_path / "full.h5"), results(tmp_path / "resumed.h5")
    assert expected.keys() == actual.keys()
    assert {"Specific Heat", "Binder Cumulant", "Magnetization^2 slope"} <= expected.keys()
    for name, reference in expected.items():
        assert reference.count == actual[name].count == 74
        for field in ("mean", "error", "covariance", "batch_sums", "batch_counts"):
            np.testing.assert_array_equal(getattr(reference, field), getattr(actual[name], field),
                                          err_msg=name + "/" + field)
    loaded = {item.props["observable"]: item
              for item in pyalps.loadMeasurements([str(tmp_path / "resumed.h5")])[0]}
    np.testing.assert_array_equal([item.mean for item in loaded["Energy"].y], actual["Energy"].mean)
    assert len(loaded["Bond-type Energy"].y) == actual["Bond-type Energy"].mean.size
    with h5py.File(tmp_path / "resumed.h5") as archive:
        for chain in range(2):
            assert archive[f"simulation/realizations/0/clones/{chain}/measurements"][()] == 37


@pytest.mark.parametrize("broken_path", ["measurements/Energy", "checkpoint/physical_moments/0"] )
def test_spinmc_later_chain_failure_preserves_outputs(executable, tmp_path, broken_path):
    stopped = run_file(tmp_path, "stopped", execution={"max_sweeps": 14},
                       output={"checkpoint": "broken.h5"})
    execute(executable, stopped)
    with h5py.File(tmp_path / "broken.h5", "a") as archive:
        del archive["simulation/realizations/0/clones/1/" + broken_path]
    resumed = run_file(tmp_path, "resumed", input={"checkpoint": "broken.h5"},
                       output={"checkpoint": "existing-checkpoint.h5"})
    (tmp_path / "resumed.h5").write_bytes(b"existing analysis")
    (tmp_path / "existing-checkpoint.h5").write_bytes(b"existing checkpoint")
    before = {path.name: path.read_bytes() for path in tmp_path.iterdir()}
    invoked = subprocess.run([executable, str(resumed)], capture_output=True, text=True, timeout=20)
    assert invoked.returncode != 0
    assert {path.name: path.read_bytes() for path in tmp_path.iterdir()} == before


def test_spinmc_resume_can_change_stopping_and_diagnostics(executable, tmp_path):
    full = run_file(tmp_path, "full")
    stopped = run_file(tmp_path, "stopped", execution={"max_sweeps": 14},
                       output={"checkpoint": "partial.h5"})
    execute(executable, [full, stopped])
    resumed = run_file(tmp_path, "resumed", input={"checkpoint": "partial.h5"},
                       execution={"error_variable": "Energy", "error_limit": 1e-30,
                                  "print_sweeps": 2})
    invoked = subprocess.run([executable, str(resumed)], capture_output=True, text=True, timeout=20)
    assert invoked.returncode == 0, invoked.stderr
    assert invoked.stdout.strip()  # Current diagnostic interval survives base load.
    expected, actual = results(tmp_path / "full.h5"), results(tmp_path / "resumed.h5")
    for name in expected:
        np.testing.assert_array_equal(expected[name].batch_sums, actual[name].batch_sums)
        np.testing.assert_array_equal(expected[name].batch_counts, actual[name].batch_counts)

    stop_early = run_file(tmp_path, "early", input={"checkpoint": "partial.h5"},
                          execution={"error_variable": "Energy", "error_limit": 1e9})
    execute(executable, stop_early)
    assert results(tmp_path / "early.h5")["Energy"].count < 74


def test_spinmc_pools_independent_unequal_partial_bins(executable, tmp_path):
    (tmp_path / "mixed.xml").write_text('''<LATTICES>
<GRAPH name="mixed" vertices="5" dimension="1">
<EDGE source="1" target="2" type="0"/><EDGE source="2" target="3" type="2"/>
<EDGE source="3" target="4" type="0"/><EDGE source="4" target="5" type="2"/>
</GRAPH></LATTICES>''', encoding="utf-8")
    parameters = {"THERMALIZATION": 0, "SWEEPS": 37, "GRAPH": "mixed", "J2": [-.3]}
    graph = {"lattice_library": "mixed.xml"}
    joint = run_file(tmp_path, "joint", parameters=parameters,
                     execution={"max_sweeps": 13}, input=graph)
    first = run_file(tmp_path, "first", parameters=parameters,
                     execution={"chains": 1, "max_sweeps": 13}, input=graph)
    second = run_file(tmp_path, "second", parameters=parameters,
                      execution={"chains": 1, "seed": 2874, "disorder_seed": 2873,
                                 "max_sweeps": 13}, input=graph)
    execute(executable, [joint, first, second])
    pooled, left, right = (results(tmp_path / (name + ".h5")) for name in ("joint", "first", "second"))
    assert pooled["Bond-type Energy"].mean.size > 1
    for name in ("Energy", "Energy^2", "Magnetization^2", "Bond-type Energy"):
        assert pooled[name].count == 26
        np.testing.assert_array_equal(pooled[name].batch_counts,
                                      np.concatenate((left[name].batch_counts, right[name].batch_counts)))
        np.testing.assert_array_equal(pooled[name].batch_sums,
                                      np.concatenate((left[name].batch_sums, right[name].batch_sums)))


@pytest.mark.parametrize("parameters", [{"MODEL": "Ising", "UPDATE": "local"}, {"MODEL": "Ising", "UPDATE": "cluster"}, {"MODEL": "O(4)", "UPDATE": "cluster"}, {"MODEL": "Potts", "q": 3, "UPDATE": "cluster"}])
def test_spinmc_extend_completed_run(executable, tmp_path, parameters):
    full = run_file(tmp_path, "full", parameters=parameters)
    short = run_file(tmp_path, "short", parameters={**parameters, "SWEEPS": 8},
                     output={"checkpoint": "continuation.h5"})
    execute(executable, [full, short])
    resumed = run_file(tmp_path, "resumed", parameters=parameters,
                       input={"checkpoint": "continuation.h5"},
                       output={"checkpoint": "extended.h5"})
    execute(executable, resumed)
    expected, actual = results(tmp_path / "full.h5"), results(tmp_path / "resumed.h5")
    assert expected.keys() == actual.keys()
    for name in expected:
        for field in ("mean", "error", "batch_sums", "batch_counts"):
            np.testing.assert_array_equal(getattr(expected[name], field), getattr(actual[name], field))
    invalid = run_file(tmp_path, "too-short", parameters={**parameters, "SWEEPS": 36},
                       input={"checkpoint": "extended.h5"},
                       output={"checkpoint": "invalid.h5"})
    before = (tmp_path / "extended.h5").read_bytes()
    rejected = subprocess.run([executable, str(invalid)], text=True, capture_output=True, timeout=20)
    assert rejected.returncode != 0
    assert (tmp_path / "extended.h5").read_bytes() == before
    assert not (tmp_path / "too-short.h5").exists()


def test_native_diagnostics_retain_chronology_and_exact_hierarchy(executable, tmp_path):
    full = run_file(tmp_path, "diagnostics", execution={"bins": 64})
    stopped = run_file(tmp_path, "stopped", execution={"bins": 64, "max_sweeps": 14},
                       output={"checkpoint": "partial.h5"})
    execute(executable, [full, stopped])
    resumed = run_file(tmp_path, "resumed", execution={"bins": 64},
                       input={"checkpoint": "partial.h5"})
    execute(executable, resumed)
    for chain in range(2):
        base = f"/simulation/realizations/0/clones/{chain}"
        with hdf5.archive(tmp_path / "diagnostics.h5") as archive:
            series = alea.BatchAccumulator.read(archive, base + "/series/Energy")
            reference = alea.AutocorrelationResult.read(archive, base + "/autocorrelation/Energy")
        with hdf5.archive(tmp_path / "resumed.h5") as archive:
            restarted = alea.AutocorrelationResult.read(archive, base + "/autocorrelation/Energy")
        bins = series.result()
        order = np.argsort(series.batch_offsets)
        order = order[bins.batch_counts[order] > 0]
        np.testing.assert_array_equal(series.batch_offsets[order], np.arange(37))
        np.testing.assert_array_equal(bins.batch_counts[order], 1)
        oracle = alea.AutocorrelationAccumulator()
        for sample in bins.batch_sums[order]:
            oracle << sample
        expected = oracle.result()
        for i in range(reference.levels):
            for field in ("mean", "variance", "error"):
                np.testing.assert_array_equal(getattr(reference.level(i), field),
                                              getattr(restarted.level(i), field))
                np.testing.assert_array_equal(getattr(reference.level(i), field),
                                              getattr(expected.level(i), field))
    curves = pyalps.loadBinningAnalysis([str(tmp_path / "diagnostics.h5")], "Energy")[0]
    assert len(curves) == 1
    assert curves[0].native_result.levels == reference.levels
    assert curves[0].y.shape == (reference.levels,)
    assert len(pyalps.loadBinningAnalysis([str(tmp_path / "diagnostics.h5")], "Energy",
        respath="/simulation/realizations/0/clones/1/autocorrelation")[0]) == 1


@pytest.mark.parametrize("parameters,unavailable", [
    ({"SWEEPS": 1}, {"Specific Heat", "Connected Susceptibility", "Binder Cumulant U2",
                     "Binder Cumulant", "Binder Cumulant slope", "Magnetization^2 slope", "Magnetization^4 slope"}),
    ({"SWEEPS": 17, "S": 0.}, {"Binder Cumulant U2", "Binder Cumulant", "Binder Cumulant slope"}),
])
def test_undefined_estimates_have_explicit_reasons(executable, tmp_path, parameters, unavailable):
    run = run_file(tmp_path, "undefined", parameters=parameters, execution={"chains": 1})
    execute(executable, run)
    with h5py.File(tmp_path / "undefined.h5") as archive:
        reasons = {pyalps.hdf5_name_decode(name): value.asstr()[()]
                   for name, value in archive["simulation/unavailable"].items()}
        assert reasons.keys() == unavailable
        assert all(reasons.values())
        assert not unavailable.intersection(pyalps.hdf5_name_decode(name)
                    for name in archive["simulation/results"])
        if parameters.get("S") == 0.:
            assert archive["simulation/results/Specific Heat/mean/value"][0] == 0.
