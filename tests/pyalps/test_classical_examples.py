"""Classical native examples: bond Hamiltonians, replica physics and restart."""
import itertools
import subprocess

import h5py
import numpy as np
import pyalps
import pytest
from pyalps.run_io import write_run_file
from conftest import tutorials_build
from test_native_mpi import compare, invoke, launcher


@pytest.fixture(params=['ising', 'heisenberg'])
def model(request):
    return request.param


def run(directory, model, name, *, exchange=True, parameters=None, budget=0,
        checkpoint=None, rng='mt19937', chains=2, processes=1, launcher=None, parallel='chains'):
    executable = str(tutorials_build() / f'00-examples/parapack/{model}/{model}')
    p = dict(ALGORITHM=model + ('; exchange' if exchange else ''), GRAPH='test',
             J=.7, J0=.8, J2=-.4, SWEEPS=127, THERMALIZATION=17)
    if exchange:
        p['INVERSE_TEMPERATURE_SET'] = [0., .3, .9]
    else:
        p['BETA'] = .3
    p.update(parameters or {})
    graph = directory/'graph.xml'
    if not graph.exists():
        edges = ([(1, 2, 0), (2, 3, 2), (3, 4, 0), (4, 1, 2), (2, 2, 0)]
                 if model == 'ising' else [(1, 2, 2)])
        graph.write_text('<LATTICES><GRAPH name="test" vertices="' + ('4' if model == 'ising' else '2')
                         + '" dimension="1">' + ''.join(f'<EDGE source="{a}" target="{b}" type="{t}"/>'
                         for a, b, t in edges) + '</GRAPH></LATTICES>')
    path = write_run_file(directory/(name+'.toml'), parameters=p,
        input=dict(lattice_library=str(graph), **({'checkpoint': str(checkpoint)} if checkpoint else {})),
        execution=dict(seed=71831, bins=16, rng=rng, chains=chains, max_sweeps=budget, parallel=parallel),
        output=dict(results=name+'.h5', checkpoint=name+'.checkpoint.h5'))
    invoke(launcher, executable, path, processes=processes)
    return directory/(name+'.h5'), directory/(name+'.checkpoint.h5')


def expected(model, beta):
    if model == 'ising':
        spins = np.array(list(itertools.product((-1., 1.), repeat=4)))
        energy = -.8*(spins[:, 0]*spins[:, 1]+spins[:, 2]*spins[:, 3]+1)
        energy += .4*(spins[:, 1]*spins[:, 2]+spins[:, 3]*spins[:, 0])
        mag = spins.mean(axis=1)
        values = dict(Energy=energy, **{'Energy^2': energy**2, 'Magnetization': mag,
                      'Magnetization^2': mag**2, 'Magnetization^4': mag**4})
        weights = np.exp(-beta*energy)
        result = {name: np.average(value, weights=weights) for name, value in values.items()}
        sites = 4
    else:
        # Two unit vectors: cos(theta) is uniform on [-1,1] before Boltzmann
        # weighting. Isotropy relates z-component moments to magnitude moments.
        c, w = np.polynomial.legendre.leggauss(80)
        energy = .4*c
        weights = w*np.exp(-beta*energy)
        average = lambda value: np.average(value, weights=weights)
        m2, m4 = average((1+c)/2), average((1+c)**2/4)
        result = {'Energy': average(energy), 'Energy^2': average(energy**2),
                  'Magnetization Z': 0., 'Magnetization^2': m2, 'Magnetization^4': m4,
                  'Magnetization Z^2': m2/3, 'Magnetization Z^4': m4/5,
                  'Binder Ratio of Magnetization Z': (m2/3)**2/(m4/5)}
        sites = 2
    result['Specific Heat'] = beta**2*(result['Energy^2']-result['Energy']**2)/sites
    result['Energy Density'] = result['Energy']/sites
    result['Binder Ratio of Magnetization'] = result['Magnetization^2']**2/result['Magnetization^4']
    return result


@pytest.mark.parametrize('exchange', [False, True])
def test_bond_hamiltonian_and_component_moments(model, tmp_path, exchange):
    filename, _ = run(tmp_path, model, 'physics', exchange=exchange, parameters=dict(SWEEPS=60000))
    groups = pyalps.loadMeasurements([str(filename)])
    assert len(groups) == (3 if exchange else 1)
    for beta, group in zip((0., .3, .9) if exchange else (.3,), groups):
        results = {d.props['observable']: d.native_result for d in group}
        for name, target in expected(model, beta).items():
            result = results[name]
            assert result.count == 120000
            np.testing.assert_allclose(result.mean, [target], atol=max(.008, 7*result.error[0]), err_msg=name)
        assert group[0].props['BETA'] == beta
        if beta == 0:
            assert 'T' not in group[0].props
            np.testing.assert_array_equal(results['Specific Heat'].mean, [0.])
        if exchange:
            assert results['EXMC: Inverse Temperature'].count == 120000


@pytest.mark.parametrize('rng', ['mt19937', 'lagged_fibonacci607'])
@pytest.mark.parametrize('mode', ['temperature', 'ladder', 'rate', 'population', 'no-exchange'])
def test_feedback_and_partial_batch_continuation(model, tmp_path, rng, mode):
    p = {} if mode == 'temperature' else dict(RANDOM_EXCHANGE=True, EXCHANGE_INTERVAL=3)
    if mode in ('rate', 'population'):
        p.update(OPTIMIZE_TEMPERATURE=True, OPTIMIZATION_TYPE=mode,
                 INITIAL_BLOCK_SWEEPS=17, OPTIMIZATION_ITERATIONS=1)
        if mode == 'population':
            p['INVERSE_TEMPERATURE_SET'] = [.1, .3, .9]
    elif mode == 'no-exchange':
        p['NO_EXCHANGE'] = True
    options = dict(exchange=mode != 'temperature', rng=rng, chains=1)
    full, full_state = run(tmp_path, model, 'full', parameters=p, **options)
    for budget in (7, 41, 107):
        _, state = run(tmp_path, model, f'part{budget}', parameters=p, budget=budget, **options)
        result, restored = run(tmp_path, model, f'restored{budget}', parameters=p, checkpoint=state, **options)
        compare(full, result); compare(full_state, restored)
    _, extended = run(tmp_path, model, 'extended', parameters=dict(p, SWEEPS=191), checkpoint=full_state, **options)
    _, reference = run(tmp_path, model, 'reference', parameters=dict(p, SWEEPS=191), **options)
    compare(reference, extended)


@pytest.mark.parametrize('fault', ['spin', 'topology', 'count', 'diagnostics', 'rng'])
def test_corrupt_checkpoint_keeps_existing_outputs(model, tmp_path, fault):
    _, state = run(tmp_path, model, 'part', budget=53)
    result, resumed = run(tmp_path, model, 'resume', checkpoint=state)
    before = result.read_bytes(), resumed.read_bytes()
    with h5py.File(state, 'r+') as ar:
        clone = ar['simulation/realizations/0/clones/1']
        if fault == 'spin': clone['walkers/0/spins'][0, :] = 0
        elif fault == 'topology': clone['walkers/0/topology'][0] = 999
        elif fault == 'count': clone['exchange/production'][()] += 1
        elif fault == 'diagnostics': del clone['stages/0/exchange/EXMC: Inverse Temperature']
        else: del clone['walkers/0/rng']
    executable = tutorials_build()/f'00-examples/parapack/{model}/{model}'
    rejected = subprocess.run([str(executable), str(tmp_path/'resume.toml')], capture_output=True, timeout=60)
    assert rejected.returncode != 0
    assert before == (result.read_bytes(), resumed.read_bytes())


def test_mpi_independent_ladders_preserve_evidence(model, tmp_path, launcher):
    full, full_state = run(tmp_path, model, 'full', chains=3)
    _, partial = run(tmp_path, model, 'part', chains=3, budget=41, processes=2, launcher=launcher)
    result, restored = run(tmp_path, model, 'resume', chains=3, checkpoint=partial, processes=3, launcher=launcher)
    compare(full, result); compare(full_state, restored)


@pytest.mark.parametrize('fault', ['temperature', 'ladder-on-fixed', 'population-zero',
                                  'uniform-temperature-zero', 'wrong-model', 'seed', 'bins'])
def test_validate_rejects_invalid_models_without_writing(model, tmp_path, fault):
    executable = tutorials_build()/f'00-examples/parapack/{model}/{model}'
    p = dict(ALGORITHM=model+'; exchange', LATTICE='chain lattice', L=4,
             INVERSE_TEMPERATURE_SET=[0., .3, .9], SWEEPS=31)
    execution = dict(seed=17, bins=8)
    if fault == 'temperature':
        p = dict(ALGORITHM=model, LATTICE='chain lattice', L=4, T=0., SWEEPS=31)
    elif fault == 'ladder-on-fixed':
        p.update(ALGORITHM=model, BETA=.3)
    elif fault == 'population-zero':
        p.update(OPTIMIZE_TEMPERATURE=True, OPTIMIZATION_TYPE='population')
    elif fault == 'uniform-temperature-zero':
        del p['INVERSE_TEMPERATURE_SET']
        p.update(NUM_REPLICAS=3, BETA_MIN=0., BETA_MAX=1., TEMPERATURE_DISTRIBUTION_TYPE=2)
    elif fault == 'wrong-model':
        p['ALGORITHM'] = 'heisenberg; exchange' if model == 'ising' else 'ising; exchange'
    elif fault == 'seed':
        execution['seed'] = 2147483646
    else:
        execution['bins'] = 3
    destination = tmp_path/'keep.h5'
    destination.write_bytes(b'existing scientific results')
    path = write_run_file(tmp_path/'invalid.toml', parameters=p, execution=execution,
                          output=dict(results='keep.h5'))
    process = subprocess.run([str(executable), '--validate', str(path)], capture_output=True, timeout=60)
    assert process.returncode != 0
    assert destination.read_bytes() == b'existing scientific results'


@pytest.mark.parametrize('mode', ['fixed', 'rate', 'population'])
@pytest.mark.parametrize('rng', ['mt19937', 'lagged_fibonacci607'])
def test_mpi_distributed_replicas_preserve_chronology(model, tmp_path, launcher, mode, rng):
    p = dict(RANDOM_EXCHANGE=True, EXCHANGE_INTERVAL=3)
    if mode != 'fixed':
        p.update(OPTIMIZE_TEMPERATURE=True, OPTIMIZATION_TYPE=mode,
                 INITIAL_BLOCK_SWEEPS=17, OPTIMIZATION_ITERATIONS=1,
                 INVERSE_TEMPERATURE_SET=[.1, .3, .9])
    options = dict(chains=2 if mode == 'fixed' else 1, parameters=p, rng=rng)
    full, full_state = run(tmp_path, model, 'full', **options)
    _, partial = run(tmp_path, model, 'part', budget=41, parallel='replicas',
                     processes=2, launcher=launcher, **options)
    result, restored = run(tmp_path, model, 'resume', checkpoint=partial,
                          parallel='replicas', processes=4, launcher=launcher, **options)
    compare(full, result); compare(full_state, restored)
    # Execution layout is not scientific state: resume a serial checkpoint
    # across MPI teams, and the MPI checkpoint back on one process.
    _, serial_partial = run(tmp_path, model, 'serial-part', budget=107, **options)
    result, restored = run(tmp_path, model, 'mpi-resume', checkpoint=serial_partial,
                          parallel='replicas', processes=3, launcher=launcher, **options)
    compare(full, result); compare(full_state, restored)
    result, restored = run(tmp_path, model, 'serial-resume', checkpoint=partial, **options)
    compare(full, result); compare(full_state, restored)


def test_mpi_fixed_temperature_rejects_replica_layout(model, tmp_path, launcher):
    executable = str(tutorials_build()/f'00-examples/parapack/{model}/{model}')
    p = dict(ALGORITHM=model, LATTICE='chain lattice', L=7, BETA=.3, SWEEPS=31)
    output = tmp_path/'keep.h5'
    output.write_bytes(b'Existing scientific output')
    path = write_run_file(tmp_path/'bad.toml', parameters=p,
                          execution=dict(parallel='replicas'), output=dict(results=output.name))
    failed = invoke(launcher, executable, path, processes=2, success=False)
    assert 'parallel = replicas requires an exchange algorithm' in failed.stdout + failed.stderr
    assert output.read_bytes() == b'Existing scientific output'


def test_mpi_replica_layout_requires_consensus(model, tmp_path, launcher):
    from test_spatial_ising_example import mpmd_failure
    executable = str(tutorials_build()/f'00-examples/parapack/{model}/{model}')
    p = dict(ALGORITHM=model+'; exchange', LATTICE='chain lattice', L=7,
             INVERSE_TEMPERATURE_SET=[0., .3, .9], SWEEPS=31)
    paths = [write_run_file(tmp_path/(mode+'.toml'), parameters=p,
                           execution=dict(parallel=mode), output=dict(results='keep.h5'))
             for mode in ('chains', 'replicas')]
    output = tmp_path/'keep.h5'
    output.write_bytes(b'Existing scientific output')
    mpmd_failure(executable, launcher, [[path] for path in paths])
    assert output.read_bytes() == b'Existing scientific output'


def test_mpi_multiple_runs_keep_separate_execution_layouts(tmp_path, launcher):
    import os
    import shlex
    executable = str(tutorials_build()/'00-examples/parapack/ising/ising')
    paths = []
    for layout in ('chains', 'replicas'):
        name = 'reference-'+layout
        run(tmp_path, 'ising', name, parallel=layout)
        path = tmp_path/(layout+'.toml')
        path.write_text((tmp_path/(name+'.toml')).read_text().replace(name, layout))
        paths.append(str(path))
    result = subprocess.run([launcher, *shlex.split(os.environ.get('ALPS_MPIEXEC_ARGS', '')),
                             '-n', '2', executable, *paths], capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stdout+result.stderr
    for layout in ('chains', 'replicas'):
        for suffix in ('.h5', '.checkpoint.h5'):
            compare(tmp_path/('reference-'+layout+suffix), tmp_path/(layout+suffix))
