# Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
"""Native Wang-Landau restart, independent-chain evaluation and exact physics."""
from pathlib import Path
import shutil
import subprocess

import h5py
import numpy as np
import pytest
from pyalps import tools
from pyalps.run_io import execute, write_run_file
from test_native_mpi import compare
from conftest import alps_program


@pytest.fixture
def executable():
    return alps_program("qwl")


def run_file(directory, name, *, parameters=None, execution=None, input=None):
    p = dict(LATTICE="chain lattice", MODEL="spin", L=4, J=1., local_S=.5,
             CUTOFF=16, NUMBER_OF_WANG_LANDAU_STEPS=4, SWEEPS=4000)
    p.update(parameters or {})
    e = dict(seed=137, bins=8)
    e.update(execution or {})
    return write_run_file(directory / (name + ".toml"), parameters=p, execution=e,
        input=input, output=dict(results=name + ".out.h5", checkpoint=name + ".checkpoint.h5"), overwrite=True)


@pytest.mark.parametrize("window,zhou,coupling,combinatorics,rng", [
    (3, True, 1., True, "lagged_fibonacci607"), (0, False, -1., False, "mt19937")])
def test_exact_restart_and_extension(executable, tmp_path, rng, window, zhou, coupling, combinatorics):
    p = dict(EXPANSION_ORDER_MINIMUM=window, USE_ZHOU_BHATT_METHOD=zhou, BLOCK_SWEEPS=200,
             J=coupling, INCLUDE_COMBINATORICS_FACTORS=combinatorics, NUMBER_OF_WANG_LANDAU_STEPS=8)
    e = dict(rng=rng, chains=2)
    full = run_file(tmp_path, "full", parameters=p, execution=e)
    partial = run_file(tmp_path, "partial", parameters=p, execution=dict(e, max_sweeps=137))
    execute(executable, [full, partial])
    resumed = run_file(tmp_path, "resumed", parameters=p, execution=e,
                       input=dict(checkpoint="partial.checkpoint.h5"))
    execute(executable, resumed)
    compare(tmp_path / "full.out.h5", tmp_path / "resumed.out.h5")
    compare(tmp_path / "full.checkpoint.h5", tmp_path / "resumed.checkpoint.h5")
    # Completed-target extension is shared by these modes; the shifted window
    # exercises it with nonzero expansion-order offsets and the alternate RNG.
    if not window:
        return
    longer = run_file(tmp_path, "longer", parameters=dict(p, SWEEPS=7000), execution=e)
    extended = run_file(tmp_path, "extended", parameters=dict(p, SWEEPS=7000), execution=e,
                        input=dict(checkpoint="full.checkpoint.h5"))
    execute(executable, [longer, extended])
    compare(tmp_path / "longer.out.h5", tmp_path / "extended.out.h5")
    compare(tmp_path / "longer.checkpoint.h5", tmp_path / "extended.checkpoint.h5")


def test_restart_during_multicanonical_production(executable, tmp_path):
    rng = "lagged_fibonacci607"
    execute(executable, run_file(tmp_path, "full", execution=dict(rng=rng)))
    state = "simulation/realizations/0/clones/0/checkpoint/"
    with h5py.File(tmp_path / "full.checkpoint.h5") as ar:
        refinement = int(ar[state + "sweeps"][()] - ar[state + "production_sweeps"][()])
    execute(executable, run_file(tmp_path, "partial", execution=dict(rng=rng, max_sweeps=refinement + 137)))
    with h5py.File(tmp_path / "partial.checkpoint.h5") as ar:
        assert ar[state + "doing_multicanonical"][()]
        assert ar[state + "production_sweeps"][()] == 137
        assert not ar[state + "all_done"][()]
    execute(executable, run_file(tmp_path, "resumed", execution=dict(rng=rng),
                                input=dict(checkpoint="partial.checkpoint.h5")))
    compare(tmp_path / "full.out.h5", tmp_path / "resumed.out.h5")
    compare(tmp_path / "full.checkpoint.h5", tmp_path / "resumed.checkpoint.h5")


def test_bad_checkpoint_preserves_outputs(executable, tmp_path):
    execute(executable, run_file(tmp_path, "partial", execution=dict(max_sweeps=13)))
    checkpoint = tmp_path / "partial.checkpoint.h5"
    original = checkpoint.read_bytes()
    for fault in ['g/values', 'state', 'missing_measurement']:
        checkpoint.write_bytes(original)
        with h5py.File(tmp_path / "partial.checkpoint.h5", "a") as ar:
            root = ar["simulation/realizations/0/clones/0"]
            if fault == "missing_measurement":
                del root["measurements/Time Up"]
            elif fault == "g/values":
                root["checkpoint/g/values"][0] = np.nan
            else:
                root["checkpoint/" + fault][...] = 999999
        run = run_file(tmp_path, "bad", input=dict(checkpoint="partial.checkpoint.h5"))
        (tmp_path / "bad.out.h5").write_bytes(b"existing results")
        (tmp_path / "bad.checkpoint.h5").write_bytes(b"existing checkpoint")
        before = {p.name: p.read_bytes() for p in tmp_path.iterdir()}
        result = subprocess.run([executable, str(run)], capture_output=True, timeout=20)
        assert result.returncode != 0
        assert {p.name: p.read_bytes() for p in tmp_path.iterdir()} == before


@pytest.mark.parametrize("coupling,combinatorics", [(-1., False)])
def test_thermodynamics_matches_exact_diagonalization(executable, tmp_path, coupling, combinatorics):
    run = run_file(tmp_path, "physics", parameters=dict(J=coupling, SWEEPS=40000,
        NUMBER_OF_WANG_LANDAU_STEPS=8, CUTOFF=24, INCLUDE_COMBINATORICS_FACTORS=combinatorics),
        execution=dict(chains=4))
    execute(executable, run)
    evaluator = str(Path(executable).with_name("qwl_evaluate"))
    result_path = tmp_path / "physics results.analysis"
    shutil.copyfile(tmp_path / "physics.out.h5", result_path)
    data = tools.evaluateQWL(result_path, appname=evaluator,
                            T_MIN=1., T_MAX=1., DELTA_T=.1)
    curves = {dataset.props["ylabel"]: dataset for group in data for dataset in group}
    # Independent dense Hamiltonian for the four-site periodic spin-1/2 chain.
    h = np.zeros((16, 16))
    magnetization = np.array([i.bit_count() - 2 for i in range(16)])
    for state in range(16):
        for site in range(4):
            neighbor = (site + 1) % 4
            same = ((state >> site) & 1) == ((state >> neighbor) & 1)
            h[state, state] += coupling * (.25 if same else -.25)
            if not same:
                h[state ^ (1 << site) ^ (1 << neighbor), state] += coupling / 2
    energies, vectors = np.linalg.eigh(h)
    weights = np.exp(-energies)
    z = weights.sum()
    energy = np.dot(energies, weights) / z
    magnetic = np.dot(np.sum(vectors**2 * magnetization[:, None]**2, axis=0), weights) / z / 4
    expected = {"Energy Density": energy / 4, "Free Energy Density": -np.log(z) / 4,
        "Entropy Density": (energy + np.log(z)) / 4,
        "Specific Heat per Site": (np.dot(energies**2, weights) / z - energy**2) / 4,
        "Uniform Structure Factor per Site": magnetic, "Uniform Susceptibility per Site": magnetic}
    for name, value in expected.items():
        assert abs(float(curves[name].y[0]) - value) < .035, (name, curves[name].y, value)


def test_validation_and_incomplete_evaluation(executable, tmp_path):
    evaluator = str(Path(executable).with_name("qwl_evaluate"))
    for parameters in (dict(CUTOFF=0), dict(EXPANSION_ORDER_MINIMUM=17), dict(local_S=1.),
                       dict(h=.2), dict(Jz=.7), dict(J=0.), dict(INITIAL_MODIFICATION_FACTOR=1.)):
        run = run_file(tmp_path, "invalid", parameters=parameters)
        result = subprocess.run([executable, "--validate", str(run)], capture_output=True, timeout=20)
        assert result.returncode != 0, parameters
        assert not (tmp_path / "invalid.out.h5").exists()
    execute(executable, run_file(tmp_path, "partial", execution=dict(max_sweeps=1)))
    result = subprocess.run([evaluator, str(tmp_path / "partial.out.h5")], capture_output=True, timeout=20)
    assert result.returncode != 0
    assert not list(tmp_path.glob("*.plot.*.xml"))


def test_histogram_completion_without_magnetic_measurements(executable, tmp_path):
    run = write_run_file(tmp_path / "automatic.toml", parameters=dict(
        LATTICE="triangular lattice", MODEL="spin", L=3, J=-1., local_S=.5,
        CUTOFF=16, NUMBER_OF_WANG_LANDAU_STEPS=3, MEASURE_MAGNETIC_PROPERTIES=False),
        execution=dict(seed=137), output=dict(results="automatic.h5"))
    execute(executable, run)
    with h5py.File(tmp_path / "automatic.h5") as ar:
        assert ar["simulation/realizations/0/clones/0/complete"][()]
        assert "Uniform Structure Factor Coefficients" not in ar["simulation/results"]
    evaluator = str(Path(executable).with_name("qwl_evaluate"))
    curves = tools.evaluateQWL(tmp_path / "automatic.h5", appname=evaluator,
                              T_MIN=1., T_MAX=1.)[0]
    assert {curve.props["ylabel"] for curve in curves} == {
        "Energy Density", "Free Energy Density", "Entropy Density", "Specific Heat per Site"}
