# Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
"""Native loop: physical checkpoints and signed estimator alignment."""
import os
from pathlib import Path
import subprocess

import h5py
import numpy as np
import pytest
from pyalps.run_io import write_run_file
from test_native_mpi import compare, invoke, launcher


@pytest.fixture
def executable():
    value = os.environ.get('ALPS_LOOP_EXECUTABLE')
    if not value:
        pytest.skip('set ALPS_LOOP_EXECUTABLE')
    return value


def run(executable, directory, name, parameters, *, rng='mt19937', budget=0,
        checkpoint=None, chains=2, bins=16):
    p = dict(ALGORITHM='loop', LATTICE='chain lattice', MODEL='spin', local_S=.5,
             L=4, J=1., T=1., THERMALIZATION=100, SWEEPS=1000)
    p.update(parameters)
    config = write_run_file(directory / (name + '.toml'), parameters=p,
        execution=dict(seed=137, bins=bins, chains=chains, rng=rng, max_sweeps=budget),
        input=dict(checkpoint=checkpoint) if checkpoint else None,
        output=dict(results=name + '.h5', checkpoint=name + '.checkpoint.h5'))
    result = subprocess.run([executable, str(config)], capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, result.stdout + result.stderr
    return directory / (name + '.checkpoint.h5')


@pytest.mark.parametrize('algorithm', ['loop', 'loop; sse'])
@pytest.mark.parametrize('parameters', [{}, {'local_S': 1.}, {'L': 3},
                                      {'L': 3, 'DISABLE_IMPROVED_ESTIMATOR': True}])
@pytest.mark.parametrize('rng', ['mt19937', 'lagged_fibonacci607'])
def test_exact_continuation(executable, tmp_path, algorithm, parameters, rng):
    p = dict(parameters, ALGORITHM=algorithm)
    run(executable, tmp_path, 'full', p, rng=rng)
    for phase, budget in [('thermal', 31), ('production', 247)]:
        run(executable, tmp_path, phase, p, rng=rng, budget=budget)
        run(executable, tmp_path, phase + '-resumed', p, rng=rng,
            checkpoint=phase + '.checkpoint.h5')
        for suffix in ('.h5', '.checkpoint.h5'):
            compare(tmp_path / ('full' + suffix), tmp_path / (phase + '-resumed' + suffix))
    run(executable, tmp_path, 'longer', dict(p, SWEEPS=1400), rng=rng)
    run(executable, tmp_path, 'extended', dict(p, SWEEPS=1400), rng=rng,
        checkpoint='full.checkpoint.h5')
    compare(tmp_path / 'longer.checkpoint.h5', tmp_path / 'extended.checkpoint.h5')


@pytest.mark.parametrize('algorithm', ['loop', 'loop; sse'])
def test_normal_estimator_sign_matches_physical_state(executable, tmp_path, algorithm):
    # Odd periodic Heisenberg chain: each off-diagonal bond vertex contributes
    # a minus sign. A nine-sweep CT trajectory contains a sign-changing flip.
    checkpoint = run(executable, tmp_path, 'signed', dict(ALGORITHM=algorithm, L=3,
        THERMALIZATION=0, SWEEPS=9, DISABLE_IMPROVED_ESTIMATOR=True), chains=32, bins=32)
    with h5py.File(checkpoint) as ar:
        for clone in ar['simulation/realizations/0/clones'].values():
            operators = clone['checkpoint/operators'][()]
            sign = (-1.) ** np.sum(operators[:, 0] & 1)
            magnetization = np.sum(.5 - clone['checkpoint/spins'][()])
            batch = clone['measurements/Magnetization^2/batch']
            occupied = np.flatnonzero(batch['count'][()])
            latest = occupied[np.argmax(batch['offset'][()][occupied])]
            assert batch['count'][latest] == 1
            np.testing.assert_array_equal(batch['sum'][latest], [sign * magnetization**2, sign])


def spin_expectations(sites, beta=1., field=0.):
    from itertools import product
    states = list(product((-.5, .5), repeat=sites))
    indices = {state: i for i, state in enumerate(states)}
    h = np.zeros((len(states), len(states)))
    for i, state in enumerate(states):
        for left in range(sites):
            right = (left + 1) % sites
            h[i, i] += state[left] * state[right]
            if state[left] != state[right]:
                target = list(state)
                target[left], target[right] = target[right], target[left]
                h[i, indices[tuple(target)]] += .5
    h -= np.diag([field * sum(state) for state in states])
    energies, vectors = np.linalg.eigh(h)
    weights = np.exp(-beta * energies)
    weights /= weights.sum()
    magnetization = np.array([sum(state) for state in states])
    expected = {'Energy': energies @ weights}
    for power in (2, 4):
        expected['Magnetization^' + str(power)] = (magnetization**power @ (vectors**2)) @ weights
    expected['Susceptibility'] = beta * expected['Magnetization^2'] / sites
    expected['Specific Heat'] = beta**2 * ((energies**2) @ weights - expected['Energy']**2) / sites
    expected['Binder Ratio of Magnetization'] = expected['Magnetization^2']**2 / expected['Magnetization^4']
    return expected


@pytest.mark.parametrize('algorithm', ['loop', 'loop; sse'])
@pytest.mark.parametrize('improved', [False, True])
@pytest.mark.parametrize('sites', [3, 4])
def test_matches_exact_diagonalization(executable, tmp_path, algorithm, improved, sites):
    expected = spin_expectations(sites)
    run(executable, tmp_path, 'physics', dict(ALGORITHM=algorithm, L=sites,
        THERMALIZATION=2000, SWEEPS=50000, DISABLE_IMPROVED_ESTIMATOR=not improved),
        chains=3, bins=64)
    with h5py.File(tmp_path / 'physics.h5') as ar:
        for name, exact in expected.items():
            result = ar['simulation/results/' + name]
            mean, error = result['mean/value'][0], result['mean/error'][0]
            assert abs(mean - exact) < max(.015, 5 * error), (name, mean, error, exact)


@pytest.mark.parametrize('fault', ['chain', 'spins', 'model', 'measurements'])
def test_failed_load_preserves_outputs(executable, tmp_path, fault):
    checkpoint = run(executable, tmp_path, 'partial', {}, budget=247)
    with h5py.File(checkpoint, 'a') as ar:
        clone = ar['simulation/realizations/0/clones/1']
        if fault == 'measurements':
            del clone['measurements/Energy']
        else:
            clone['checkpoint/' + fault][...] = 255
    config = write_run_file(tmp_path / 'bad.toml', parameters=dict(
        ALGORITHM='loop', LATTICE='chain lattice', MODEL='spin', local_S=.5,
        L=4, J=1., T=1., THERMALIZATION=100, SWEEPS=1000),
        execution=dict(seed=137, bins=16, chains=2),
        input=dict(checkpoint=checkpoint.name),
        output=dict(results='keep.h5', checkpoint='keep.checkpoint.h5'))
    for name in ('keep.h5', 'keep.checkpoint.h5'):
        (tmp_path / name).write_bytes(b'original')
    before = {p.name: p.read_bytes() for p in tmp_path.iterdir()}
    result = subprocess.run([executable, str(config)], capture_output=True, timeout=60)
    assert result.returncode != 0
    assert before == {p.name: p.read_bytes() for p in tmp_path.iterdir()}


REPLICA_MODES = [
    {}, {'RANDOM_EXCHANGE': True, 'EXCHANGE_INTERVAL': 3}, {'NO_EXCHANGE': True},
    {'OPTIMIZE_TEMPERATURE': True, 'OPTIMIZATION_TYPE': 'rate', 'EXCHANGE_INTERVAL': 3},
    {'OPTIMIZE_TEMPERATURE': True, 'OPTIMIZATION_TYPE': 'population'},
]


@pytest.mark.parametrize('algorithm,rng', [('loop; exchange', 'mt19937'),
                                         ('loop; sse; exchange', 'lagged_fibonacci607')])
@pytest.mark.parametrize('mode', REPLICA_MODES)
def test_replica_continuation_and_loader(executable, tmp_path, algorithm, rng, mode):
    import pyalps
    p = dict(ALGORITHM=algorithm, NUM_REPLICAS=3, T_MIN=.8, T_MAX=1.2,
             INITIAL_BLOCK_SWEEPS=100, OPTIMIZATION_ITERATIONS=1, **mode)
    run(executable, tmp_path, 'full', p, chains=1, bins=8, rng=rng)
    phases = [('feedback', 31), ('production', 757)]
    if mode.get('OPTIMIZE_TEMPERATURE'):
        phases.append(('optimization', 151))
    for phase, budget in phases:
        run(executable, tmp_path, phase, p, chains=1, bins=8, rng=rng, budget=budget)
        run(executable, tmp_path, phase + '-resumed', p, chains=1, bins=8, rng=rng,
            checkpoint=phase + '.checkpoint.h5')
        for suffix in ('.h5', '.checkpoint.h5'):
            compare(tmp_path / ('full' + suffix), tmp_path / (phase + '-resumed' + suffix))
    groups = pyalps.loadMeasurements([str(tmp_path / 'full.h5')], ['Energy'])
    assert len(groups) == 3
    for index, group in enumerate(groups):
        assert len(group) == 1
        assert group[0].props['replica'] == index
        assert group[0].native_result.count == 1000
    temperatures = [group[0].props['T'] for group in groups]
    assert temperatures[0] == pytest.approx(1.2)
    assert temperatures[-1] == pytest.approx(.8)
    assert temperatures[0] > temperatures[1] > temperatures[2]
    diagnostics = pyalps.loadBinningAnalysis([str(tmp_path / 'full.h5')], ['Energy'])
    assert len(diagnostics) == 3
    assert [group[0].props['T'] for group in diagnostics] == temperatures
    with h5py.File(tmp_path / 'full.h5') as ar:
        assert 'results' not in ar['simulation']


@pytest.mark.parametrize('algorithm', ['loop; exchange', 'loop; sse; exchange'])
@pytest.mark.parametrize('rng', ['mt19937', 'lagged_fibonacci607'])
def test_checkpoint_offdiagonal_validity(executable, tmp_path, algorithm, rng):
    # Graph decorations are chosen before a cluster flip. For an offdiagonal
    # XXZ vertex their diagonal predicate need not match the saved worldline.
    # The SSE kernel has no longitudinal-field implementation.
    p = dict(ALGORITHM=algorithm, L=3, h=0. if 'sse' in algorithm else .4, DISABLE_IMPROVED_ESTIMATOR=True,
             NO_EXCHANGE=True, NUM_REPLICAS=3, T_MIN=.8, T_MAX=1.2,
             SWEEPS=127, THERMALIZATION=17)
    full = run(executable, tmp_path, 'full', p, bins=8, rng=rng)
    partial = run(executable, tmp_path, 'partial', p, bins=8, rng=rng, budget=31)
    resumed = run(executable, tmp_path, 'resumed', p, bins=8, rng=rng,
                  checkpoint=partial.name)
    compare(full, resumed)
    compare(tmp_path / 'full.h5', tmp_path / 'resumed.h5')

    with h5py.File(partial, 'a') as ar:
        candidates = [replica['checkpoint']
                      for clone in ar['simulation/realizations/0/clones'].values()
                      for replica in clone['replicas'].values()]
        for state in candidates:
            operators = state['operators'][()]
            offdiagonal = np.flatnonzero((operators[:, 0] & 1) & (operators[:, 1] == 1))
            if not offdiagonal.size:
                continue
            # Two flips on the same parallel pair are periodic but have zero
            # XXZ matrix element. Periodicity alone must not accept this state.
            bad = operators[[offdiagonal[0], offdiagonal[0]]]
            state['spins'][...] = 0
            attributes = dict(state['operators'].attrs)
            del state['operators']
            dataset = state.create_dataset('operators', data=bad)
            for key, value in attributes.items(): dataset.attrs[key] = value
            if 'times' in state:
                attributes = dict(state['times'].attrs)
                del state['times']
                dataset = state.create_dataset('times', data=[.25, .75])
                for key, value in attributes.items(): dataset.attrs[key] = value
            break
        else:
            pytest.fail('Checkpoint fixture needs a bond offdiagonal operator')
    outputs = tmp_path / 'resumed.h5', resumed
    before = [path.read_bytes() for path in outputs]
    result = subprocess.run([executable, str(tmp_path / 'resumed.toml')],
                            text=True, capture_output=True, timeout=60)
    assert result.returncode != 0
    assert 'operator incompatible with worldline spins' in result.stderr
    assert before == [path.read_bytes() for path in outputs]


@pytest.mark.parametrize('algorithm,sites,field,mode', [
    ('loop; exchange', 3, 0., {}),
    ('loop; sse; exchange', 3, 0., {'RANDOM_EXCHANGE': True}),
    ('loop; exchange', 4, .7, {'OPTIMIZE_TEMPERATURE': True, 'OPTIMIZATION_TYPE': 'rate'}),
    ('loop; exchange', 4, 0., {'OPTIMIZE_TEMPERATURE': True, 'OPTIMIZATION_TYPE': 'population'}),
])
def test_replica_physics(executable, tmp_path, algorithm, sites, field, mode):
    import pyalps
    run(executable, tmp_path, 'physics', dict(ALGORITHM=algorithm, L=sites, h=field,
        NUM_REPLICAS=3, T_MIN=.6, T_MAX=1.7, THERMALIZATION=2000, SWEEPS=30000,
        INITIAL_BLOCK_SWEEPS=1000, OPTIMIZATION_ITERATIONS=1, **mode), chains=1, bins=64)
    groups = pyalps.loadMeasurements([str(tmp_path / 'physics.h5')])
    assert len(groups) == 3
    for group in groups:
        values = {d.props['observable']: d.native_result for d in group}
        expected = spin_expectations(sites, 1 / group[0].props['T'], field)
        for name, exact in expected.items():
            value = values[name]
            assert abs(value.mean[0] - exact) < max(.025, 5 * value.error[0]), (name, value.mean, value.error, exact)


@pytest.mark.parametrize('algorithm,rng,optimize', [
    ('loop; exchange', 'mt19937', False), ('loop; sse; exchange', 'lagged_fibonacci607', True),
])
def test_replica_mpi_restart(executable, launcher, tmp_path, algorithm, rng, optimize):
    p = dict(ALGORITHM=algorithm, LATTICE='chain lattice', MODEL='spin', local_S=.5,
        L=4, J=1., TEMPERATURE_SET=[.8, 1., 1.2], THERMALIZATION=100, SWEEPS=500,
        OPTIMIZE_TEMPERATURE=optimize, OPTIMIZATION_TYPE='population',
        INITIAL_BLOCK_SWEEPS=100, OPTIMIZATION_ITERATIONS=1)
    def config(name, budget=0, checkpoint=None):
        return write_run_file(tmp_path / (name + '.toml'), parameters=p,
            execution=dict(seed=137, bins=8, chains=1 if optimize else 3, rng=rng, max_sweeps=budget),
            input=dict(checkpoint=checkpoint) if checkpoint else None,
            output=dict(results=name + '.h5', checkpoint=name + '.checkpoint.h5'))
    invoke(launcher, executable, config('serial'))
    invoke(launcher, executable, config('mpi'), processes=2)
    invoke(launcher, executable, config('partial', 31), processes=2)
    invoke(launcher, executable, config('resumed', checkpoint='partial.checkpoint.h5'), processes=3)
    for name in ('mpi', 'resumed'):
        for suffix in ('.h5', '.checkpoint.h5'):
            compare(tmp_path / ('serial' + suffix), tmp_path / (name + suffix))


def test_replica_common_hamiltonian(executable, tmp_path):
    import pyalps
    checkpoint = run(executable, tmp_path, 'common', dict(ALGORITHM='loop; exchange', J='T',
        T=1., NUM_REPLICAS=3, T_MIN=.6, T_MAX=1.7, THERMALIZATION=2000, SWEEPS=30000),
        chains=1, bins=64)
    with h5py.File(checkpoint) as ar:
        walkers = ar['simulation/realizations/0/clones/0/replicas']
        for index in ('1', '2'):
            np.testing.assert_array_equal(walkers['0/checkpoint/model'][()], walkers[index + '/checkpoint/model'][()])
    for group in pyalps.loadMeasurements([str(tmp_path / 'common.h5')], ['Energy']):
        value = group[0].native_result
        exact = spin_expectations(4, 1 / group[0].props['T'])['Energy']
        assert abs(value.mean[0] - exact) < max(.025, 5 * value.error[0])


@pytest.mark.parametrize('fault', ['beta', 'walkers', 'weights', 'production'])
def test_replica_corrupt_checkpoint(executable, tmp_path, fault):
    p = dict(ALGORITHM='loop; exchange', LATTICE='chain lattice', MODEL='spin', local_S=.5,
        L=4, J=1., T=1., THERMALIZATION=100, SWEEPS=1000, NUM_REPLICAS=3, T_MIN=.8, T_MAX=1.2)
    checkpoint = run(executable, tmp_path, 'partial', p, budget=247, chains=1, bins=8)
    with h5py.File(checkpoint, 'a') as ar:
        value = ar['simulation/realizations/0/clones/0/exchange/' + fault]
        value[...] = np.nan if fault == 'weights' else 0
    config = write_run_file(tmp_path / 'bad.toml', parameters=p,
        execution=dict(seed=137, bins=8, chains=1), input=dict(checkpoint=checkpoint.name),
        output=dict(results='keep.h5', checkpoint='keep.checkpoint.h5'))
    for name in ('keep.h5', 'keep.checkpoint.h5'):
        (tmp_path / name).write_bytes(b'original')
    before = {p.name: p.read_bytes() for p in tmp_path.iterdir()}
    assert subprocess.run([executable, str(config)], capture_output=True, timeout=60).returncode != 0
    assert before == {p.name: p.read_bytes() for p in tmp_path.iterdir()}
