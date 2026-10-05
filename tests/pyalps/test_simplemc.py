# Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
"""The installed simplemc driver must preserve native statistics and restart."""
import os
from pathlib import Path
import shutil
import subprocess

import h5py
import numpy as np
import pytest

import pyalps
from pyalps import alea, hdf5, tools
from pyalps.run_io import execute, write_run_file


@pytest.fixture
def executable():
    explicit = os.environ.get("ALPS_SIMPLEMC_EXECUTABLE")
    if explicit:
        return str(Path(explicit).resolve(strict=True))
    try:
        tools.check_existence("simplemc")
    except RuntimeError:
        pytest.skip("native simplemc executable is not installed")
    return shutil.which("simplemc")


def run_file(directory, name, model="ising", *, execution=None, input=None, output=None, parameters=None):
    p = dict(ALGORITHM=model, LATTICE="chain lattice", L=5, T=1.8,
             H=.17, J=.8, SWEEPS=37, THERMALIZATION=3)
    p.update(parameters or {})
    e = dict(seed=2873, chains=2, bins=8)
    e.update(execution or {})
    o = dict(results=name + ".h5")
    o.update(output or {})
    return write_run_file(directory / (name + ".toml"), parameters=p,
                          execution=e, input=input, output=o)


def read_results(filename):
    with hdf5.archive(filename) as archive:
        return {pyalps.hdf5_name_decode(name): alea.BatchResult.read(
                    archive, "/simulation/results/" + name)
                for name in archive.list_children("/simulation/results")}


@pytest.mark.parametrize("model", ["ising", "xy", "heisenberg"])
@pytest.mark.parametrize("rng", ["mt19937", "lagged_fibonacci607"])
def test_simplemc_restart_matches_complete_native_batches(executable, tmp_path, model, rng):
    full = run_file(tmp_path, "full", model, execution={"rng": rng}, output={"checkpoint": "full-checkpoint.h5"})
    stopped = run_file(tmp_path, "stopped", model, execution={"max_sweeps": 14, "rng": rng},
                       output={"checkpoint": "stopped-checkpoint.h5"})
    execute(executable, [full, stopped])
    resumed = run_file(tmp_path, "resumed", model, execution={"rng": rng}, input={"checkpoint": "stopped-checkpoint.h5"})
    execute(executable, resumed)
    expected, actual = read_results(tmp_path / "full.h5"), read_results(tmp_path / "resumed.h5")
    assert expected.keys() == actual.keys()
    assert "Specific Heat" in expected and "Binder Ratio of Magnetization" in expected
    for name in expected:
        assert expected[name].count == actual[name].count == 74
        for field in ("mean", "error", "covariance", "batch_sums", "batch_counts"):
            np.testing.assert_array_equal(getattr(expected[name], field), getattr(actual[name], field),
                                          err_msg=model + ": " + name + "/" + field)
    measured = {entry.props["observable"]: entry
                for entry in pyalps.loadMeasurements([str(tmp_path / "resumed.h5")])[0]}
    np.testing.assert_array_equal([value.mean for value in measured["Energy"].y], actual["Energy"].mean)
    np.testing.assert_array_equal([value.error for value in measured["Energy"].y], actual["Energy"].error)
    with h5py.File(tmp_path / "resumed.h5") as archive:
        for chain in range(2):
            base = f"simulation/realizations/0/clones/{chain}"
            assert archive[base + "/completed_sweeps"][()] == 40
            assert archive[base + "/measurements"][()] == 37


@pytest.mark.parametrize("thermalization,cut", [(0, 2), (3, 2)])
def test_simplemc_zero_and_interrupted_thermalization(executable, tmp_path, thermalization, cut):
    full = run_file(tmp_path, "full", parameters={"THERMALIZATION": thermalization})
    stopped = run_file(tmp_path, "stopped", parameters={"THERMALIZATION": thermalization},
                       execution={"max_sweeps": cut}, output={"checkpoint": "stopped-checkpoint.h5"})
    execute(executable, [full, stopped])
    with h5py.File(tmp_path / "stopped.h5") as archive:
        assert archive["simulation/realizations/0/clones/0/measurements"][()] == max(0, cut - thermalization)
    resumed = run_file(tmp_path, "resumed", parameters={"THERMALIZATION": thermalization},
                       input={"checkpoint": "stopped-checkpoint.h5"})
    execute(executable, resumed)
    for name, expected in read_results(tmp_path / "full.h5").items():
        actual = read_results(tmp_path / "resumed.h5")[name]
        np.testing.assert_array_equal(actual.batch_sums, expected.batch_sums)
        np.testing.assert_array_equal(actual.batch_counts, expected.batch_counts)


def test_simplemc_time_limit_keeps_resumable_state(executable, tmp_path):
    stopped = run_file(tmp_path, "stopped", execution={"time_limit": 1e-12},
                       output={"checkpoint": "stopped-checkpoint.h5"})
    execute(executable, stopped)
    with h5py.File(tmp_path / "stopped.h5") as archive:
        assert archive["simulation/realizations/0/clones/0/completed_sweeps"][()] == 0
    full = run_file(tmp_path, "full")
    resumed = run_file(tmp_path, "resumed", input={"checkpoint": "stopped-checkpoint.h5"})
    execute(executable, [full, resumed])
    np.testing.assert_array_equal(read_results(tmp_path / "full.h5")["Energy"].batch_sums,
                                  read_results(tmp_path / "resumed.h5")["Energy"].batch_sums)


def test_simplemc_malformed_later_chain_preserves_output(executable, tmp_path):
    stopped = run_file(tmp_path, "stopped", execution={"max_sweeps": 14},
                       output={"checkpoint": "broken.h5"})
    execute(executable, stopped)
    with h5py.File(tmp_path / "broken.h5", "a") as archive:
        del archive["simulation/realizations/0/clones/1/measurements/Energy"]
    resumed = run_file(tmp_path, "resumed", input={"checkpoint": "broken.h5"})
    destination = tmp_path / "resumed.h5"
    destination.write_bytes(b"existing scientific output")
    before = {p.name: p.read_bytes() for p in tmp_path.iterdir()}
    result = subprocess.run([executable, str(resumed)], text=True, capture_output=True)
    assert result.returncode != 0
    assert {p.name: p.read_bytes() for p in tmp_path.iterdir()} == before


def test_simplemc_custom_graph_bindings_couplings_and_vtk(executable, tmp_path):
    library = tmp_path / "custom.xml"
    library.write_text('''<LATTICES>
<LATTICE name="line" dimension="1"><BASIS><VECTOR>spacing</VECTOR></BASIS></LATTICE>
<UNITCELL name="cell" dimension="1"><VERTEX/>
<EDGE type="2"><SOURCE vertex="1" offset="0"/><TARGET vertex="1" offset="1"/></EDGE></UNITCELL>
<LATTICEGRAPH name="custom"><FINITELATTICE><LATTICE ref="line"/>
<EXTENT dimension="1" size="site_count"/><BOUNDARY type="periodic"/></FINITELATTICE>
<UNITCELL ref="cell"/></LATTICEGRAPH></LATTICES>''', encoding="utf-8")
    file = write_run_file(tmp_path / "custom.toml",
        parameters={"ALGORITHM": "ising", "LATTICE": "custom", "site_count": 4, "spacing": 1.7,
                    "J": 7., "J2": .3, "T": 2., "SWEEPS": 8, "THERMALIZATION": 0},
        input={"lattice_library": "custom.xml"},
        output={"results": "custom.h5", "snapshot_prefix": "spins"},
        execution={"seed": 9, "chains": 1, "bins": 4, "snapshot_interval": 8})
    execute(executable, file)
    result = read_results(tmp_path / "custom.h5")["Energy"]
    assert result.count == 8 and abs(result.mean[0]) <= 1.2 + 1e-12
    text = (tmp_path / "spins.clone1.8.vtk").read_text()
    assert "# vtk DataFile Version" in text and "POINTS 4 " in text and "SCALARS spins" in text


@pytest.mark.parametrize("parameters", [{"ALGORITHM": model} for model in ("ising", "xy", "heisenberg")])
def test_simplemc_extend_completed_run(executable, tmp_path, parameters):
    full = run_file(tmp_path, "full", parameters=parameters)
    short = run_file(tmp_path, "short", parameters={**parameters, "SWEEPS": 8},
                     output={"checkpoint": "continuation.h5"})
    execute(executable, [full, short])
    resumed = run_file(tmp_path, "resumed", parameters=parameters,
                       input={"checkpoint": "continuation.h5"},
                       output={"checkpoint": "extended.h5"})
    execute(executable, resumed)
    expected, actual = read_results(tmp_path / "full.h5"), read_results(tmp_path / "resumed.h5")
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
