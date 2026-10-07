"""Native replacements for the three scheduler Ising examples."""
import itertools
import subprocess

import h5py
import numpy as np
import pyalps
import pytest
from pyalps.run_io import write_run_file
from test_native_mpi import compare
from conftest import tutorials_build


@pytest.fixture(params=['ising1', 'ising2', 'ising3'])
def executable(request):
    return tutorials_build()/'00-examples/scheduler'/request.param


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


def test_exact_continuation(executable, tmp_path):
    # The variants share the driver; spread warmup and measurement stops and both RNGs over them.
    budget, rng = {'ising1': (9, 'mt19937'), 'ising2': (43, 'lagged_fibonacci607'),
                   'ising3': (43, 'mt19937')}[executable.name]
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
    executable = tutorials_build()/'00-examples/scheduler/ising2'
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


def test_native_binder_analysis(executable, tmp_path):
    from pyalps import alea, hdf5
    # The two evaluators differ only in the cumulant's factor; one variant each suffices.
    command, factor = ('evaluate2', 3) if executable.name == 'ising3' else ('evaluate', 1)
    source, _ = run(executable, tmp_path, 'source')
    original = source.read_bytes()
    path = tmp_path/'analysis.toml'
    output = tmp_path/'analysis.h5'
    write_run_file(path, input={'results': str(source)}, output={'results': str(output)})
    evaluator = executable.parent/command
    process = subprocess.run([str(evaluator), '--validate', str(path)], capture_output=True)
    if executable.name == 'ising2':
        assert process.returncode != 0  # This variant has no m²/m⁴ joint samples.
        assert not output.exists() and source.read_bytes() == original
        return
    assert process.returncode == 0, process.stderr
    assert not output.exists()
    process = subprocess.run([str(evaluator), str(path)], capture_output=True)
    assert process.returncode == 0, process.stderr
    assert source.read_bytes() == original
    with hdf5.archive(str(source)) as archive:
        joint = alea.read_result(archive, '/simulation/joint')
    data = {d.props['observable']: d.native_result for d in pyalps.loadMeasurements([str(output)])[0]}
    result = data['Binder cumulant of Magnetization']
    keep = joint.batch_counts > 0
    weights, sums = joint.batch_counts[keep].astype(float), joint.batch_sums[keep]
    total = sums.sum(axis=0)
    count = weights.sum()
    mean = total/count
    leave = (total-sums)/(count-weights[:, None])
    pseudo = factor*(count*mean[3]/mean[2]**2-(count-weights)*leave[:, 3]/leave[:, 2]**2)
    expected = pseudo.sum()/count
    variance = np.sum(weights*(pseudo/weights-expected)**2)/(count-np.sum(weights**2)/count)
    np.testing.assert_allclose(result.batch_sums[keep, 0], pseudo, rtol=1e-10)
    np.testing.assert_allclose(result.mean, [expected], rtol=1e-12)
    np.testing.assert_allclose(result.error, [np.sqrt(variance*np.sum(weights**2)/count**2)], rtol=1e-10)
    np.testing.assert_allclose(data['Correlations'].mean, joint.mean[4:], atol=1e-14)
    with hdf5.archive(str(output)) as archive:
        saved = alea.read_result(archive, '/simulation/joint')
    np.testing.assert_array_equal(saved.batch_counts, joint.batch_counts)
    np.testing.assert_array_equal(saved.batch_sums, joint.batch_sums)
    # Archive alias protection happens before publication.
    path.write_text('[input]\nresults = '+repr(str(source))+'\n[output]\nresults = '+repr(str(source))+'\n')
    assert subprocess.run([str(evaluator), str(path)], capture_output=True).returncode != 0
    assert source.read_bytes() == original
    # An interrupted warmup has no analysable batches; keep existing output.
    empty, _ = run(executable, tmp_path, 'warmup', budget=9)
    before = output.read_bytes()
    write_run_file(path, input={'results': str(empty)}, output={'results': str(output)}, overwrite=True)
    assert subprocess.run([str(evaluator), str(path)], capture_output=True).returncode != 0
    assert output.read_bytes() == before


def test_hot_two_site_chain_is_ergodic(executable, tmp_path):
    # At this temperature exp(-delta/T) rounds to one. Exactly two forced
    # flips per recorded sweep would preserve spin parity and freeze energy.
    path = tmp_path/'hot.toml'
    parameters = dict(L=2, T=1.e308, SWEEPS=20000, THERMALIZATION=0)
    if executable.name == 'ising2':
        parameters['LATTICE'] = 'chain lattice'
    write_run_file(path, parameters=parameters, execution=dict(seed=19, chains=1, bins=32),
                   output=dict(results='hot.h5'))
    process = subprocess.run([str(executable), str(path)], capture_output=True)
    assert process.returncode == 0, process.stderr
    data = {d.props['observable']: d.native_result for d in pyalps.loadMeasurements([str(tmp_path/'hot.h5')])[0]}
    np.testing.assert_allclose(data['Energy'].mean, [0.], atol=.025)
    if executable.name != 'ising2':
        np.testing.assert_allclose(data['Magnetization^2'].mean, [.5], atol=.025)
        np.testing.assert_allclose(data['Magnetization^4'].mean, [.5], atol=.025)
