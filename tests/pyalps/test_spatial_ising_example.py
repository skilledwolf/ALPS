"""Spatial MPI shares one physical chain and preserves serial evidence."""
import itertools
import os
from pathlib import Path
import shlex
import subprocess

import h5py
import numpy as np
import pyalps
import pytest
from pyalps import alea, hdf5
from pyalps.run_io import write_run_file
from test_native_mpi import compare, invoke, launcher
from test_single_ising_example import run
from conftest import tutorials_build, unavailable


@pytest.fixture
def executable():
    return str(tutorials_build()/'00-examples/parapack/multiple/ising_multiple')


@pytest.fixture
def mpi_executable():
    # Never infer this from the serial build: launching a serial program under
    # mpiexec could create several competing writers and conceal missing MPI.
    value = os.environ.get('ALPS_SPATIAL_ISING_MPI_EXECUTABLE')
    if not value:
        if os.environ.get('ALPS_MPIEXEC'):
            unavailable('set ALPS_SPATIAL_ISING_MPI_EXECUTABLE to the MPI-enabled example')
        pytest.skip('set ALPS_SPATIAL_ISING_MPI_EXECUTABLE to an MPI-enabled build')
    return str(Path(value).resolve(strict=True))


@pytest.mark.parametrize('scan', [False, True])
def test_serial_arbitrary_graph_fallback(executable, tmp_path, scan):
    single = str(tutorials_build()/'00-examples/parapack/single/ising_single')
    reference, reference_state = run(single, tmp_path, 'single', scan=scan, lattice='square lattice', length=4)
    actual, actual_state = run(executable, tmp_path, 'spatial', scan=scan, lattice='square lattice', length=4)
    compare(reference, actual)
    compare(reference_state, actual_state)


def test_large_ring_validation_preserves_supported_sizes(mpi_executable, launcher, tmp_path):
    path = write_run_file(tmp_path/'large.toml', parameters=dict(LATTICE='chain lattice', L=65536,
                          T=2., SWEEPS=2, THERMALIZATION=0), output=dict(results='large.h5'))
    command = [launcher, *shlex.split(os.environ.get('ALPS_MPIEXEC_ARGS', '')), '-n', '2',
               mpi_executable, '--validate', str(path)]
    validated = subprocess.run(command, capture_output=True, text=True, timeout=60)
    assert validated.returncode == 0, validated.stdout+validated.stderr
    assert not (tmp_path/'large.h5').exists()


@pytest.mark.parametrize('rng', ['mt19937', 'lagged_fibonacci607'])
@pytest.mark.parametrize('scan', [False, True])
def test_mpi_partition_independent_evidence(mpi_executable, launcher, tmp_path, rng, scan):
    reference, reference_state = run(mpi_executable, tmp_path, 'one', rng=rng, scan=scan,
                                     length=7, coupling=.7)
    for processes in (2, 3):
        actual, actual_state = run(mpi_executable, tmp_path, 'rank'+str(processes), rng=rng,
                                  scan=scan, length=7, coupling=.7, processes=processes, launcher=launcher)
        compare(reference, actual)
        compare(reference_state, actual_state)
    with h5py.File(reference_state) as archive:
        clones = archive['simulation/realizations/0/clones']
        assert len(clones) == 2
        for clone in clones.values():
            spins = clone['checkpoint/spins'][()]
            assert spins.shape == (7,) and np.all(np.isin(spins, (-1, 1)))
    groups = pyalps.loadMeasurements([str(reference)])
    assert len(groups) == (3 if scan else 1)
    for group in groups:
        values = {d.props['observable']: d.native_result for d in group}
        assert values['Energy'].count == 234
        assert abs(values['Energy'].mean[0]) <= .7*7
        assert abs(values['Magnetization'].mean[0]) <= 7
        assert 0 <= values['Magnetization^2'].mean[0] <= 49


@pytest.mark.parametrize('rng', ['mt19937', 'lagged_fibonacci607'])
def test_mpi_threaded_ring_preserves_evidence(mpi_executable, launcher, openmp_examples, tmp_path, monkeypatch, rng):
    # Three 512-site blocks give each color 256 sites, exercising the
    # threaded kernel rather than its small-lattice serial shortcut.
    monkeypatch.setenv('OMP_NUM_THREADS', '1')
    reference, reference_state = run(mpi_executable, tmp_path, 'threads1', rng=rng,
                                     length=1536, coupling=.7, sweeps=37,
                                     processes=3, launcher=launcher)
    for processes, threads in ((3, 2), (3, 4), (1, 4)):
        monkeypatch.setenv('OMP_NUM_THREADS', str(threads))
        actual, actual_state = run(mpi_executable, tmp_path, f'threads{threads}ranks{processes}', rng=rng,
                                  length=1536, coupling=.7, sweeps=37,
                                  processes=processes, launcher=launcher)
        compare(reference, actual)
        compare(reference_state, actual_state)


@pytest.mark.parametrize('rng', ['mt19937', 'lagged_fibonacci607'])
@pytest.mark.parametrize('scan,budget', [(False, 5), (False, 38), (True, 9),
                                        (True, 71), (True, 142), (True, 301)])
def test_mpi_cross_rank_continuation(mpi_executable, launcher, tmp_path, rng, scan, budget):
    reference, reference_state = run(mpi_executable, tmp_path, 'full', rng=rng, scan=scan,
                                     length=7, coupling=.7)
    _, partial = run(mpi_executable, tmp_path, 'part', rng=rng, scan=scan, budget=budget,
                     length=7, coupling=.7, processes=2, launcher=launcher)
    resumed, resumed_state = run(mpi_executable, tmp_path, 'resumed', rng=rng, scan=scan,
                                 checkpoint=partial, length=7, coupling=.7, processes=3, launcher=launcher)
    compare(reference, resumed)
    compare(reference_state, resumed_state)


@pytest.mark.parametrize('rng', ['mt19937', 'lagged_fibonacci607'])
def test_mpi_completed_target_extension(mpi_executable, launcher, tmp_path, rng):
    reference, reference_state = run(mpi_executable, tmp_path, 'full', rng=rng, length=7, coupling=.7)
    _, short = run(mpi_executable, tmp_path, 'short', rng=rng, length=7, coupling=.7,
                   sweeps=23, processes=2, launcher=launcher)
    extended, extended_state = run(mpi_executable, tmp_path, 'extended', rng=rng,
                                   checkpoint=short, length=7, coupling=.7, processes=3, launcher=launcher)
    compare(reference, extended)
    compare(reference_state, extended_state)


@pytest.mark.parametrize('processes', [2, 3])
@pytest.mark.parametrize('coupling', [1., -1.])
def test_mpi_exact_ring_thermodynamics(mpi_executable, launcher, tmp_path, processes, coupling):
    filename, _ = run(mpi_executable, tmp_path, 'physics', length=7, coupling=coupling,
                      sweeps=50000, processes=processes, launcher=launcher)
    spins = np.array(list(itertools.product((-1., 1.), repeat=7)))
    energy = -coupling*np.sum(spins*np.roll(spins, 1, axis=1), axis=1)
    mag = spins.sum(axis=1)
    means = np.average(np.array([np.full(len(spins), 7.), energy, energy**2,
                                mag, mag**2, mag**4]).T, axis=0, weights=np.exp(-energy/2.))
    expected = dict(zip(('Number of Sites', 'Energy', 'Energy^2', 'Magnetization',
                         'Magnetization^2', 'Magnetization^4'), means))
    expected['Specific Heat'] = (means[2]-means[1]**2)/(4.*7)
    expected['Energy Density'] = means[1]/7
    expected['Binder Ratio of Magnetization'] = means[4]**2/means[5]
    values = {d.props['observable']: d.native_result for d in pyalps.loadMeasurements([str(filename)])[0]}
    assert values.keys() == expected.keys()
    for name, exact in expected.items():
        assert values[name].count == 100000
        np.testing.assert_allclose(values[name].mean, [exact],
                                   atol=max(.015*max(1., abs(exact)), 6*values[name].error[0]))


def protected_failure(mpi_executable, launcher, path, processes, outputs):
    before = tuple(output.read_bytes() for output in outputs)
    failure = invoke(launcher, mpi_executable, path, processes=processes, success=False)
    assert failure.stderr or failure.stdout
    assert before == tuple(output.read_bytes() for output in outputs)


def mpmd_failure(mpi_executable, launcher, arguments):
    command = [launcher, *shlex.split(os.environ.get('ALPS_MPIEXEC_ARGS', ''))]
    for rank, args in enumerate(arguments):
        if rank:
            command.append(':')
        command.extend(['-n', '1', mpi_executable, *map(str, args)])
    failure = subprocess.run(command, capture_output=True, text=True, timeout=60)
    assert failure.returncode != 0, failure.stdout+failure.stderr
    assert 'Collective ranks require identical run configurations' in failure.stdout+failure.stderr


@pytest.mark.parametrize('field,value', [('L', 7), ('chains', 2), ('rng', 'lagged_fibonacci607')])
@pytest.mark.parametrize('validate', [False, True])
def test_mpi_rank_local_valid_runs_require_consensus(mpi_executable, launcher, tmp_path, field, value, validate):
    outputs = tmp_path/'keep.h5', tmp_path/'keep.checkpoint.h5'
    for output in outputs:
        output.write_bytes(b'Existing scientific output')
    before = tuple(output.read_bytes() for output in outputs)
    parameters = dict(LATTICE='chain lattice', L=6 if field == 'L' else 7, J=.7,
                      T=2., SWEEPS=117, THERMALIZATION=11)
    execution = dict(seed=823, chains=1, bins=8, rng='mt19937')
    paths = []
    for rank in range(2):
        p, e = parameters.copy(), execution.copy()
        if rank:
            (p if field == 'L' else e)[field] = value
        path = write_run_file(tmp_path/('rank'+str(rank)+'.toml'), parameters=p, execution=e,
                              output=dict(results=outputs[0].name, checkpoint=outputs[1].name))
        valid = subprocess.run([mpi_executable, '--validate', str(path)],
                               capture_output=True, text=True, timeout=60)
        assert valid.returncode == 0, valid.stdout+valid.stderr
        paths.append(path)
    mpmd_failure(mpi_executable, launcher, [(['--validate'] if validate else [])+[path] for path in paths])
    assert before == tuple(output.read_bytes() for output in outputs)


def test_mpi_rank_local_command_mode_requires_consensus(mpi_executable, launcher, tmp_path):
    outputs = tmp_path/'keep.h5', tmp_path/'keep.checkpoint.h5'
    for output in outputs:
        output.write_bytes(b'Existing scientific output')
    before = tuple(output.read_bytes() for output in outputs)
    path = write_run_file(tmp_path/'same.toml', parameters=dict(LATTICE='chain lattice', L=7, J=.7,
                          T=2., SWEEPS=117, THERMALIZATION=11), execution=dict(seed=823, chains=1, bins=8),
                          output=dict(results=outputs[0].name, checkpoint=outputs[1].name))
    mpmd_failure(mpi_executable, launcher, [['--validate', path], [path]])
    assert before == tuple(output.read_bytes() for output in outputs)


@pytest.mark.parametrize('lattice,length', [('square lattice', 4), ('chain lattice', 5)])
def test_mpi_invalid_decomposition_preserves_outputs(mpi_executable, launcher, tmp_path, lattice, length):
    outputs = tmp_path/'keep.h5', tmp_path/'keep.checkpoint.h5'
    for output in outputs:
        output.write_bytes(b'Existing scientific output')
    path = write_run_file(tmp_path/'bad.toml', parameters=dict(LATTICE=lattice, L=length, J=1.,
                         T=2., SWEEPS=117, THERMALIZATION=11), execution=dict(seed=823, chains=2, bins=8),
                         output=dict(results=outputs[0].name, checkpoint=outputs[1].name))
    protected_failure(mpi_executable, launcher, path, 3, outputs)


@pytest.mark.parametrize('fault', ['spin', 'topology', 'count', 'rng'])
def test_mpi_corrupt_checkpoint_preserves_outputs(mpi_executable, launcher, tmp_path, fault):
    _, partial = run(mpi_executable, tmp_path, 'part', length=7, coupling=.7,
                     budget=38, processes=2, launcher=launcher)
    outputs = run(mpi_executable, tmp_path, 'resumed', length=7, coupling=.7,
                  checkpoint=partial, processes=3, launcher=launcher)
    with h5py.File(partial, 'r+') as archive:
        checkpoint = archive['simulation/realizations/0/clones/1/checkpoint']
        if fault == 'spin':
            checkpoint['spins'][0] = 0
        elif fault == 'topology':
            checkpoint['topology'][1] = 999
        elif fault == 'count':
            checkpoint['sweeps'][()] += 1
        else:
            del checkpoint['engine']
    protected_failure(mpi_executable, launcher, tmp_path/'resumed.toml', 3, outputs)


def test_mpi_collective_stopping_and_publication(mpi_executable, launcher, tmp_path):
    parameters = dict(LATTICE='chain lattice', L=7, J=.7, T=2., SWEEPS=10000000, THERMALIZATION=11)
    execution = dict(seed=823, bins=8, chains=2, time_limit=.04, checkpoint_interval=.002)
    path = write_run_file(tmp_path/'timed.toml', parameters=parameters, execution=execution,
                          output=dict(results='timed.h5', checkpoint='timed.checkpoint.h5'))
    invoke(launcher, mpi_executable, path, processes=3)
    with h5py.File(tmp_path/'timed.checkpoint.h5') as archive:
        clones = archive['simulation/realizations/0/clones']
        assert len(clones) == 2
        steps = [int(clone['checkpoint/sweeps'][()]) for clone in clones.values()]
        assert any(steps) and all(0 <= count < 10000011 for count in steps)
    with hdf5.archive(str(tmp_path/'timed.h5')) as archive:
        result = alea.read_result(archive, '/simulation/joint')
        assert result.count == sum(max(0, count-11) for count in steps)
    assert {item.name for item in tmp_path.iterdir()} == {'timed.toml', 'timed.h5', 'timed.checkpoint.h5'}
    execution.update(time_limit=0., max_sweeps=37)
    resume = write_run_file(tmp_path/'continue.toml', parameters=parameters, execution=execution,
                            input=dict(checkpoint='timed.checkpoint.h5'),
                            output=dict(results='continued.h5', checkpoint='continued.checkpoint.h5'))
    invoke(launcher, mpi_executable, resume, processes=2)
    with h5py.File(tmp_path/'continued.checkpoint.h5') as archive:
        after = [int(clone['checkpoint/sweeps'][()]) for clone in archive['simulation/realizations/0/clones'].values()]
        assert after == [count+37 for count in steps]
