"""The scheduler Wolff lesson uses native physical and statistical continuation."""
import itertools
import subprocess

import h5py
import numpy as np
import pytest
import pyalps
from pyalps.run_io import write_run_file
from test_native_mpi import compare, invoke, launcher
from conftest import tutorials_build


@pytest.fixture
def executable():
    return str(tutorials_build()/'08-alpsize/09-scheduler/wolff')


def run(executable, directory, name, *, budget=0, checkpoint=None, rng='mt19937', sweeps=101, processes=1, launcher=None):
    run_file = directory/(name+'.toml')
    write_run_file(run_file, parameters=dict(LATTICE='chain lattice', L=4, T=2.2,
                                             SWEEPS=sweeps, THERMALIZATION=11),
                   input={'checkpoint': str(checkpoint)} if checkpoint else {},
                   execution=dict(seed=823, rng=rng, chains=3, bins=8, max_sweeps=budget),
                   output=dict(results=name+'.h5', checkpoint=name+'.checkpoint.h5'))
    invoke(launcher, executable, run_file, processes=processes)
    return directory/(name+'.h5'), directory/(name+'.checkpoint.h5')


@pytest.mark.parametrize('rng', ['mt19937', 'lagged_fibonacci607'])
@pytest.mark.parametrize('budget', [5, 12, 38])
def test_exact_continuation(executable, tmp_path, rng, budget):
    full, full_state = run(executable, tmp_path, 'full', rng=rng)
    part, state = run(executable, tmp_path, 'part', rng=rng, budget=budget)
    resumed, resumed_state = run(executable, tmp_path, 'resumed', rng=rng, checkpoint=state)
    compare(full, resumed)
    compare(full_state, resumed_state)
    loaded = pyalps.loadMeasurements([str(resumed)])[0]
    assert len(loaded) == 4 and all(d.native_result.count == 303 for d in loaded)
    if budget == 5:
        with h5py.File(part) as archive:
            assert not list(archive['simulation/results'])
            assert 'Statistics' in archive['simulation/unavailable']


def test_completed_target_extension(executable, tmp_path):
    full, full_state = run(executable, tmp_path, 'full')
    _, state = run(executable, tmp_path, 'short', sweeps=19)
    extended, extended_state = run(executable, tmp_path, 'extended', checkpoint=state)
    compare(full, extended)
    compare(full_state, extended_state)


@pytest.mark.parametrize('fault', ['spin', 'neighbors', 'count', 'chain', 'rng'])
def test_invalid_checkpoint_preserves_outputs(executable, tmp_path, fault):
    _, state = run(executable, tmp_path, 'part', budget=38)
    with h5py.File(state, 'r+') as archive:
        clone = archive['simulation/realizations/0/clones/1']
        if fault == 'spin': clone['checkpoint/spins'][0] = 0
        elif fault == 'neighbors': clone['checkpoint/neighbors'][1] = 100
        elif fault == 'count': clone['checkpoint/sweeps'][()] += 1
        elif fault == 'chain': clone['checkpoint/chain_id'][()] = 0
        else:
            # Corrupt the RNG without touching physical/statistical fields.
            del clone['checkpoint/engine']
    results, checkpoint = tmp_path/'failed.h5', tmp_path/'failed.checkpoint.h5'
    results.write_bytes(b'keep results'); checkpoint.write_bytes(b'keep checkpoint')
    run_file = tmp_path/'failed.toml'
    write_run_file(run_file, parameters=dict(LATTICE='chain lattice', L=4, T=2.2, SWEEPS=101, THERMALIZATION=11),
                   input={'checkpoint': str(state)}, execution=dict(seed=823, chains=3, bins=8),
                   output=dict(results=str(results), checkpoint=str(checkpoint)))
    process = subprocess.run([executable, str(run_file)], capture_output=True)
    assert process.returncode != 0
    assert results.read_bytes() == b'keep results'
    assert checkpoint.read_bytes() == b'keep checkpoint'


def test_mpi_checkpoint_repartition(executable, launcher, tmp_path):
    full, full_state = run(executable, tmp_path, 'full')
    _, state = run(executable, tmp_path, 'part', budget=38, processes=2, launcher=launcher)
    resumed, resumed_state = run(executable, tmp_path, 'resumed', checkpoint=state)
    compare(full, resumed)
    compare(full_state, resumed_state)


def test_pooled_chain_physics(executable, tmp_path):
    result, _ = run(executable, tmp_path, 'physics', sweeps=50000)
    values = []
    for config in itertools.product((-1, 1), repeat=4):
        spins = np.array(config)
        m = spins.mean()
        values.append((-np.sum(spins*np.roll(spins, 1)), m, m*m, m**4))
    values = np.array(values)
    expected = np.average(values[:, 1:], axis=0, weights=np.exp(-values[:, 0]/2.2))
    data = {d.props['observable']: d.native_result for d in pyalps.loadMeasurements([str(result)])[0]}
    for i, name in enumerate(('Magnetization', 'Magnetization^2', 'Magnetization^4')):
        assert data[name].count == 150000
        np.testing.assert_allclose(data[name].mean, expected[i:i+1], atol=.015)
