#!/usr/bin/env python3
"""Check native measurements, C++ callbacks, and exact subclass continuation."""
from pathlib import Path
import copy
import shutil
import tempfile
import h5py
import numpy as np
import ising
import pyalps.hdf5 as hdf5
import pyalps.ngs as ngs
from pyalps.alea import BatchResult

PARAMETERS = {'L': 7, 'THERMALIZATION': 11, 'SWEEPS': 701, 'T': 2., 'SEED': 42}
NAMES = ['Correlations', 'Energy', 'Magnetization', 'Magnetization^2', 'Magnetization^4']
CLONE = '/simulation/realizations/0/clones/0'


def assert_results_equal(actual, expected):
    assert sorted(actual) == sorted(expected) == NAMES
    for name in NAMES:
        assert isinstance(actual[name], BatchResult)
        assert actual[name].count == expected[name].count
        for field in ('mean', 'error', 'covariance', 'batch_sums', 'batch_counts'):
            np.testing.assert_array_equal(getattr(actual[name], field), getattr(expected[name], field))


simulation = ising.sim(PARAMETERS)
assert isinstance(simulation, ngs.mcbase)
assert sorted(simulation.measurements) == NAMES
assert simulation.run(lambda: False)
assert simulation.sweeps == PARAMETERS['THERMALIZATION'] + PARAMETERS['SWEEPS']
results = ngs.collectResults(simulation)
assert isinstance(results, dict)
for name, result in results.items():
    assert result.count == PARAMETERS['SWEEPS']
    assert result.mean.shape == (PARAMETERS['L'] if name == 'Correlations' else 1,)
# The exact finite periodic Ising-chain partition function gives its energy.
u = np.tanh(1 / PARAMETERS['T'])
energy = -(u + u**(PARAMETERS['L'] - 1)) / (1 + u**PARAMETERS['L'])
assert abs(results['Energy'].mean[0] - energy) <= 6 * results['Energy'].error[0]

with tempfile.TemporaryDirectory() as directory:
    checkpoint = Path(directory) / 'ising.clone0.h5'
    stopped = ising.sim(PARAMETERS)
    assert not stopped.run(lambda: stopped.sweeps == 148)
    snapshot = ngs.collectResults(stopped)
    with hdf5.archive(checkpoint, 'w') as archive:
        archive.set_context('/caller')
        archive[CLONE] = stopped
        assert archive.context == '/caller'
    with hdf5.archive(checkpoint, 'r') as archive:
        assert archive['/parameters/format'] == 'alps.params.v2'
        assert sorted(archive.list_children(CLONE + '/measurements')) == NAMES
        assert archive[CLONE + '/measurements/Energy/@kind'] == 6
        level = int(archive[CLONE + '/measurements/Energy/cursor/level'])
        counts = archive[CLONE + '/measurements/Energy/batch/count']
        assert level >= 2 and np.any((counts > 0) & (counts < 2**level))
    restored = ising.sim(dict(PARAMETERS, T=1., SEED=7))
    with hdf5.archive(checkpoint, 'r') as archive:
        archive.set_context(CLONE)
        restored.load(archive)
        assert archive.context == CLONE
    assert restored.parameters == PARAMETERS
    assert restored.sweeps == stopped.sweeps
    np.testing.assert_array_equal(restored.spins, stopped.spins)
    assert_results_equal(ngs.collectResults(restored), snapshot)
    assert restored.run(lambda: False)
    np.testing.assert_array_equal(restored.spins, simulation.spins)
    assert_results_equal(ngs.collectResults(restored), results)
    assert restored.random() == simulation.random()
    # A result is a snapshot, independently owned by the caller.
    assert snapshot['Energy'].count == 137
    target = ising.sim(dict(PARAMETERS, T=1.5, SEED=31))
    assert not target.run(lambda: target.sweeps == 31)
    original = ngs.collectResults(target)
    parameters, spins, sweeps = dict(target.parameters), target.spins.copy(), target.sweeps
    handle = target.measurements['Energy']
    random = copy.deepcopy(target.random)
    random_values = [random() for _ in range(17)]
    broken = Path(directory) / 'invalid-app.h5'
    for fault in ('spins', 'sweeps', 'temperature', 'geometry', 'count'):
        shutil.copyfile(checkpoint, broken)
        with h5py.File(broken, 'a') as archive:
            if fault == 'spins':
                archive[CLONE + '/checkpoint/spins'][0] = 0
            elif fault == 'sweeps':
                archive[CLONE + '/checkpoint/sweeps'][()] = PARAMETERS['THERMALIZATION'] + PARAMETERS['SWEEPS'] + 1
            elif fault == 'count':
                archive[CLONE + '/measurements/Energy/batch/count'][0] += 1
            else:
                key = 'T' if fault == 'temperature' else 'L'
                for entry in archive['parameters/entries'].values():
                    if entry['name'].asstr()[()] == key:
                        entry['value'][()] = 0
        with hdf5.archive(broken, 'r') as archive:
            archive.set_context(CLONE)
            try:
                target.load(archive)
            except ValueError:
                pass
            else:
                raise AssertionError(f'invalid {fault} checkpoint accepted')
            assert archive.context == CLONE
        assert target.parameters == parameters and target.sweeps == sweeps
        assert target.measurements['Energy'] is handle
        np.testing.assert_array_equal(target.spins, spins)
        assert_results_equal(ngs.collectResults(target), original)
        random = copy.deepcopy(target.random)
        assert [random() for _ in range(17)] == random_values
    output = Path(directory) / 'results.h5'
    hdf5.save_checkpoint(str(output), lambda archive:
                         ngs.saveResults(results, restored.parameters, archive, '/simulation/results'))
    with hdf5.archive(output, 'r') as archive:
        assert sorted(archive.list_children('/simulation/results')) == NAMES
        assert archive['/simulation/results/Energy/@kind'] == 5

print('native Python subclass and exact checkpoint continuation: ok')
