"""Installed C++ lesson: TOML CLI, native evidence and independent physics oracle."""
import itertools
import os
from pathlib import Path
import subprocess

import numpy as np
import pytest
import pyalps
from pyalps import alea, hdf5


@pytest.fixture
def executable():
    path = os.environ.get('ALPS_CPP_ISING_EXECUTABLE')
    if not path:
        pytest.skip('Set ALPS_CPP_ISING_EXECUTABLE to the built C++ tutorial')
    return str(Path(path).resolve())


def run_file(tmp_path, parameters='', output='result.h5'):
    path = tmp_path/'run.toml'
    path.write_text('[parameters]\nL=2\nSWEEPS=100000\nTHERMALIZATION=1000\n'
                    + parameters + '\n[execution]\nseed=823\nbins=128\n'
                    + f'[output]\nresults="{output}"\n')
    return path


def test_native_evidence_and_exact_small_lattice(executable, tmp_path):
    run = run_file(tmp_path, 'BETA=0.3')
    subprocess.run([executable, str(run)], check=True, capture_output=True, text=True)
    filename = str(tmp_path/'result.h5')
    with hdf5.archive(filename) as archive:
        joint = alea.read_result(archive, '/simulation/joint')
        assert joint.count == 100000
    values = []
    for spins in itertools.product((-1, 1), repeat=4):
        spins = np.array(spins).reshape(2, 2)
        e = -np.sum(spins*(np.roll(spins, 1, 0)+np.roll(spins, 1, 1)))/4
        m = spins.mean()
        values.append([e, m, abs(m), m*m, m**4])
    values = np.array(values)
    expected = np.average(values, axis=0, weights=np.exp(-0.3*4*values[:, 0]))
    np.testing.assert_allclose(joint.mean, expected, atol=0.025)
    loaded = pyalps.loadMeasurements([filename])[0]
    assert all(d.props['L'] == 2 for d in loaded)
    measurements = {d.props['observable']: d.native_result for d in loaded}
    for i, name in enumerate(('E', 'm', '|m|', 'm^2', 'm^4')):
        result = measurements[name]
        np.testing.assert_allclose(result.batch_sums, joint.batch_sums[:, i:i+1], atol=1e-9)
        np.testing.assert_array_equal(result.batch_counts, joint.batch_counts)
        np.testing.assert_allclose(result.error, joint.error[i:i+1], rtol=1e-12)
    diagnostics = pyalps.loadBinningAnalysis([filename])[0]
    assert len(diagnostics) == 5
    assert all(d.native_result.count == 100000 for d in diagnostics)


def test_cli_preflight(executable, tmp_path):
    run = run_file(tmp_path)
    before = run.read_bytes()
    subprocess.run([executable, '--validate', str(run)], check=True, capture_output=True)
    assert run.read_bytes() == before and not (tmp_path/'result.h5').exists()
    schema = subprocess.run([executable, '--schema'], check=True, capture_output=True, text=True).stdout
    assert 'application = "cpp-ising"' in schema
    for args in (['--validate'], ['--Tmin', '1'], [str(run), str(run)]):
        assert subprocess.run([executable, *args], capture_output=True).returncode != 0
    assert not (tmp_path/'result.h5').exists()
    run.write_text(run.read_text().replace('bins=128', 'bins=3'))
    assert subprocess.run([executable, '--validate', str(run)], capture_output=True).returncode != 0
    run_file(tmp_path, output='run.toml')
    before = run.read_bytes()
    assert subprocess.run([executable, str(run)], capture_output=True).returncode != 0
    assert run.read_bytes() == before


@pytest.mark.parametrize('beta', [-0.1, 1e308])
def test_beta_boundaries(executable, tmp_path, beta):
    run = run_file(tmp_path, f'BETA={beta}')
    process = subprocess.run([executable, str(run)], capture_output=True, text=True)
    if beta < 0:
        assert process.returncode != 0 and not (tmp_path/'result.h5').exists()
    else:
        assert process.returncode == 0, process.stderr
        data = pyalps.loadMeasurements([str(tmp_path/'result.h5')], ['E'])[0][0]
        np.testing.assert_array_equal(data.native_result.mean, [-2.])
        assert np.isfinite(data.native_result.error).all()


def test_original_no_argument_scan(executable, tmp_path):
    subprocess.run([executable], cwd=tmp_path, check=True, capture_output=True)
    files = sorted(tmp_path.glob('ising.L_16beta_*.h5'))
    assert len(files) == 11
    assert (tmp_path/'ising.L_16beta_0.h5').exists()
    assert (tmp_path/'ising.L_16beta_1.h5').exists()
    datasets = pyalps.loadMeasurements([str(path) for path in files])
    assert all(len(run) == 5 for run in datasets)
    for run in datasets:
        for data in run:
            assert data.props['L'] == 16 and data.native_result.count == 5000
            assert np.isfinite(data.native_result.mean).all()
            assert np.isfinite(data.native_result.error).all()


def test_infinite_temperature_samples_all_parities(executable, tmp_path):
    path = run_file(tmp_path, 'BETA=0.0')
    subprocess.run([executable, str(path)], check=True, capture_output=True)
    with hdf5.archive(str(tmp_path/'result.h5')) as archive:
        joint = alea.read_result(archive, '/simulation/joint')
    np.testing.assert_allclose(joint.mean, [0., 0., .375, .25, .15625], atol=.02)
