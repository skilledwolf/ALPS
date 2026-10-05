"""The introductory Python lesson uses native joint statistical evidence."""
import importlib.util
import itertools
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

import pyalps
from pyalps import alea, hdf5, run_io


LESSON = Path(__file__).resolve().parents[2]/'tutorials/09-code/01-python'
spec = importlib.util.spec_from_file_location('ising_lesson', LESSON/'solution/ising.py')
ising = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ising)


def test_separate_instances_own_rng_and_samples():
    first = ising.Simulation(0.3, 2, seed=127, bins=8)
    other = ising.Simulation(0.3, 2, seed=777, bins=8)
    reference = ising.Simulation(0.3, 2, seed=127, bins=8)
    other.run(17, 91)
    first.run(13, 67)
    reference.run(13, 67)
    assert first.spins == reference.spins
    np.testing.assert_array_equal(first.samples.result().batch_sums, reference.samples.result().batch_sums)
    assert first.samples.count == 67
    assert all(acc.count == 67 for acc in first.diagnostics.values())
    with pytest.raises(ValueError, match='fresh simulation'):
        first.run(0, 2)


def test_joint_binder_error_retains_covariance(tmp_path, monkeypatch):
    sim = ising.Simulation(0., 2, bins=16)
    # Two correlated moment columns with m^4 = 3*m^2 - 2 on every sample.
    # Independent-error propagation would substantially overestimate the error.
    samples = np.tile([[0., 1., 1., 1., 1.], [0., np.sqrt(2), np.sqrt(2), 2., 4.]], (8, 1))
    sequence = iter(samples)
    monkeypatch.setattr(sim, 'step', lambda: None)
    monkeypatch.setattr(sim, 'observables', lambda: next(sequence))
    sim.run(0, len(samples))
    results, unavailable = sim.results()
    assert not unavailable
    ratio = results['Binder Ratio']
    n = len(samples)
    mean = samples.mean(axis=0)
    leave_out = (samples.sum(axis=0) - samples)/(n-1)
    estimates = leave_out[:, 4]/leave_out[:, 3]**2
    expected = n*mean[4]/mean[3]**2 - (n-1)*estimates.mean()
    error = np.sqrt((n-1)*np.mean((estimates-estimates.mean())**2))
    np.testing.assert_allclose(ratio.mean, [expected], rtol=1e-12)
    np.testing.assert_allclose(ratio.error, [error], rtol=1e-12)
    filename = tmp_path/'result.h5'
    sim.save(filename)
    loaded = pyalps.loadMeasurements([str(filename)], ['Binder Ratio'])[0][0]
    np.testing.assert_array_equal(loaded.native_result.batch_sums, ratio.batch_sums)
    with hdf5.archive(filename) as archive:
        joint = alea.read_result(archive, '/simulation/joint')
        diagnostics = alea.read_result(archive, '/simulation/realizations/0/clones/0/autocorrelation/E')
    np.testing.assert_array_equal(joint.batch_sums, sim.samples.result().batch_sums)
    assert diagnostics.count == n
    curves = pyalps.loadBinningAnalysis([str(filename)])[0]
    assert len(curves) == len(ising.NAMES)
    assert all(curve.native_result.count == n for curve in curves)


def test_undefined_binder_ratio_is_explicit():
    sim = ising.Simulation(0., 2)
    for _ in range(8):
        sim.samples << np.zeros(5)
    results, unavailable = sim.results()
    assert 'Binder Ratio' not in results
    assert 'positive <m^2>' in unavailable['Binder Ratio']


def test_small_lattice_matches_exact_enumeration():
    sim = ising.Simulation(0.3, 2, seed=823, bins=128)
    sim.run(1000, 100000)
    values = []
    for spins in itertools.product((-1, 1), repeat=4):
        sim.spins = np.array(spins).reshape(2, 2).tolist()
        values.append(sim.observables())
    values = np.array(values)
    weights = np.exp(-0.3*4*values[:, 0])
    expected = np.average(values, axis=0, weights=weights)
    np.testing.assert_allclose(sim.samples.result().mean, expected, atol=0.025)


@pytest.mark.parametrize('script', ['solution/ising.py', 'solution/ising_binder.py', 'solution/run.py', 'ising-skeleton.py'])
def test_toml_validation_keeps_inputs_and_outputs_untouched(tmp_path, script):
    run = tmp_path/'run.toml'
    run_io.write_run_file(run, ising.SCHEMA, parameters={'L': 2, 'BETA': 0.2, 'SWEEPS': 17},
                          output={'results': 'result.h5'})
    result = subprocess.run([sys.executable, str(LESSON/script), '--validate', str(run)],
                            capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stderr
    assert 'Valid Ising' in result.stdout
    assert sorted(p.name for p in tmp_path.iterdir()) == ['run.toml']


@pytest.mark.parametrize('script', ['ising.py', 'ising_binder.py'])
def test_cli_output_loads_and_collision_is_rejected(tmp_path, script):
    run = tmp_path/'run.toml'
    run_io.write_run_file(run, ising.SCHEMA,
                         parameters={'L': 2, 'BETA': 0.2, 'SWEEPS': 67, 'THERMALIZATION': 3},
                         execution={'seed': 841, 'bins': 8}, output={'results': 'result.h5'})
    executable = [sys.executable, str(LESSON/'solution'/script)]
    completed = subprocess.run(executable + [str(run)], capture_output=True, text=True, timeout=30,
                               env={**os.environ, 'MPLBACKEND': 'Agg'})
    assert completed.returncode == 0, completed.stderr
    assert 'unavailable (insufficient bins)' in completed.stdout
    data = pyalps.loadMeasurements([str(tmp_path/'result.h5')])[0]
    assert {d.props['observable'] for d in data} == set(ising.NAMES) | {'Binder Ratio'}
    assert all(d.native_result.count == 67 for d in data)
    before = run.read_bytes()
    rejected = subprocess.run(executable + [str(run), str(run)], capture_output=True, text=True, timeout=30)
    assert rejected.returncode != 0 and 'distinct' in rejected.stderr
    assert run.read_bytes() == before


def test_infinite_temperature_samples_all_parities():
    sim = ising.Simulation(0., 2, seed=823, bins=64)
    sim.run(0, 30000)
    # All sixteen spin configurations must contribute, not one parity sector.
    np.testing.assert_allclose(sim.samples.result().mean, [0., 0., .375, .25, .15625], atol=.02)
