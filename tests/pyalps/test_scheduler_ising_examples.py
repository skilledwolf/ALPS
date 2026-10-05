"""Native replacements for the three scheduler Ising examples."""
import itertools
import os
from pathlib import Path
import subprocess

import h5py
import numpy as np
import pyalps
import pytest
from pyalps.run_io import write_run_file
from test_native_mpi import compare


@pytest.fixture(params=['ising1', 'ising2', 'ising3'])
def executable(request):
    root = os.environ.get('ALPS_EXAMPLES_BUILD_DIR')
    if not root:
        pytest.skip('Set ALPS_EXAMPLES_BUILD_DIR to the built C++ examples')
    return Path(root).resolve()/'scheduler'/request.param


def run(executable, directory, name, *, sweeps=117, budget=0, checkpoint=None, rng='mt19937'):
    parameters = dict(L=4, T=2., SWEEPS=sweeps, THERMALIZATION=17)
    if executable.name == 'ising2':
        parameters['LATTICE'] = 'chain lattice'
    path = directory/(name+'.toml')
    write_run_file(path, parameters=parameters,
                   input={'checkpoint': str(checkpoint)} if checkpoint else {},
                   execution=dict(seed=19, chains=2, bins=8, max_sweeps=budget, rng=rng),
                   output=dict(results=name+'.h5', checkpoint=name+'.checkpoint.h5'))
    process = subprocess.run([str(executable), str(path)], capture_output=True, text=True)
    assert process.returncode == 0, process.stdout+process.stderr
    return directory/(name+'.h5'), directory/(name+'.checkpoint.h5')


@pytest.mark.parametrize('rng', ['mt19937', 'lagged_fibonacci607'])
@pytest.mark.parametrize('budget', [9, 43])
def test_exact_continuation(executable, tmp_path, rng, budget):
    full, full_state = run(executable, tmp_path, 'full', rng=rng)
    _, state = run(executable, tmp_path, 'part', budget=budget, rng=rng)
    resumed, resumed_state = run(executable, tmp_path, 'resumed', checkpoint=state, rng=rng)
    compare(full, resumed)
    compare(full_state, resumed_state)


def test_extension_and_invalid_physical_state(executable, tmp_path):
    full, full_state = run(executable, tmp_path, 'full')
    _, state = run(executable, tmp_path, 'short', sweeps=23)
    extended, extended_state = run(executable, tmp_path, 'extended', checkpoint=state)
    compare(full, extended)
    compare(full_state, extended_state)
    with h5py.File(state, 'r+') as archive:
        archive['simulation/realizations/0/clones/1/checkpoint/spins'][0] = 0
    before = extended.read_bytes(), extended_state.read_bytes()
    process = subprocess.run([str(executable), str(tmp_path/'extended.toml')], capture_output=True)
    assert process.returncode != 0
    assert before == (extended.read_bytes(), extended_state.read_bytes())


def test_exact_chain_physics_and_correlations(executable, tmp_path):
    path, _ = run(executable, tmp_path, 'physics', sweeps=50000)
    data = {d.props['observable']: d.native_result for d in pyalps.loadMeasurements([str(path)])[0]}
    spins = np.array(list(itertools.product((-1., 1.), repeat=4)))
    energy = -np.sum(spins*np.roll(spins, 1, axis=1), axis=1)
    magnetization = spins.mean(axis=1)
    weights = np.exp(-energy/2.)
    exact = {'Energy': np.average(energy/4., weights=weights),
             'Magnetization': np.average(magnetization, weights=weights)}
    if executable.name != 'ising2':
        exact.update({'Magnetization^2': np.average(magnetization**2, weights=weights),
                      'Magnetization^4': np.average(magnetization**4, weights=weights),
                      'Correlations': [np.average((spins*np.roll(spins, d, axis=1)).mean(axis=1),
                                                   weights=weights) for d in range(4)]})
    assert data.keys() == exact.keys()
    for name, value in exact.items():
        assert data[name].count == 100000
        np.testing.assert_allclose(data[name].mean, np.atleast_1d(value), atol=.015)
    if 'Correlations' in data:
        np.testing.assert_array_equal(data['Correlations'].mean[[0]], [1.])
        np.testing.assert_allclose(data['Correlations'].mean[1], -data['Energy'].mean[0], atol=1e-14)


def test_lattice_square_and_validation(tmp_path):
    root = os.environ.get('ALPS_EXAMPLES_BUILD_DIR')
    if not root:
        pytest.skip('Set ALPS_EXAMPLES_BUILD_DIR to the built C++ examples')
    executable = Path(root).resolve()/'scheduler/ising2'
    path = tmp_path/'square.toml'
    write_run_file(path, parameters=dict(LATTICE='square lattice', L=2, T=3.,
                                         SWEEPS=100000, THERMALIZATION=1000),
                   execution=dict(seed=74, bins=32), output=dict(results='square.h5'))
    process = subprocess.run([str(executable), '--validate', str(path)], capture_output=True)
    assert process.returncode == 0, process.stderr
    assert not (tmp_path/'square.h5').exists()
    process = subprocess.run([str(executable), str(path)], capture_output=True)
    assert process.returncode == 0, process.stderr
    spins = np.array(list(itertools.product((-1., 1.), repeat=4))).reshape(-1, 2, 2)
    energy = -(spins*(np.roll(spins, 1, axis=1)+np.roll(spins, 1, axis=2))).sum(axis=(1, 2))
    expected = np.average(energy/4., weights=np.exp(-energy/3.))
    data = {d.props['observable']: d.native_result for d in pyalps.loadMeasurements([str(tmp_path/'square.h5')])[0]}
    np.testing.assert_allclose(data['Energy'].mean, [expected], atol=.015)
    path.write_text(path.read_text().replace('3.0', '0.0'))
    before = (tmp_path/'square.h5').read_bytes()
    process = subprocess.run([str(executable), '--validate', str(path)], capture_output=True)
    assert process.returncode != 0
    assert (tmp_path/'square.h5').read_bytes() == before
