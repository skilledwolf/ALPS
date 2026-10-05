# Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
"""Production quantum kernels: exact continuation and independent physics."""
import os
from pathlib import Path
import subprocess

import h5py
import numpy as np
import pytest
from pyalps.run_io import write_run_file
from test_native_mpi import compare


@pytest.fixture
def executables():
    names = {app: os.environ.get('ALPS_' + app.upper() + '_EXECUTABLE')
             for app in ('worm', 'dirloop_sse')}
    if not all(names.values()):
        pytest.skip('set ALPS_WORM_EXECUTABLE and ALPS_DIRLOOP_SSE_EXECUTABLE')
    return names


def run(exe, directory, name, p, *, rng='mt19937', budget=0, checkpoint=None):
    parameters = dict(LATTICE='chain lattice', MODEL='spin', L=4, local_S=.5,
                      J=1., T=1., THERMALIZATION=100, SWEEPS=1000, SKIP=3)
    parameters.update(p)
    if Path(exe).name == 'dirloop_sse':
        parameters.setdefault('INITIAL_CUTOFF', 64)
    config = write_run_file(directory / (name + '.toml'), parameters=parameters,
        execution=dict(seed=137, bins=8, chains=2, rng=rng, max_sweeps=budget),
        input=dict(checkpoint=checkpoint) if checkpoint else None,
        output=dict(results=name + '.h5', checkpoint=name + '.checkpoint.h5'), overwrite=True)
    result = subprocess.run([exe, str(config)], capture_output=True, text=True, timeout=90)
    assert result.returncode == 0, result.stdout + result.stderr
    return directory / (name + '.h5')


CASES = [
    ('worm', {}),
    ('worm', dict(MODEL='boson Hubbard', Nmax=2, U=1., t=.3, mu=.5, NONLOCAL=False)),
    ('worm', dict(LATTICE='square lattice', L=2, MODEL='boson Hubbard', Nmax=2,
                  U=1., t=.3, mu=.5, CHAIN_KAPPA=True)),
    ('worm', dict(USE_1D_STIFFNESS=True)),
    ('dirloop_sse', dict(WHICH_LOOP_TYPE='minbounce')),
    ('dirloop_sse', dict(WHICH_LOOP_TYPE='heatbath', local_S=1.)),
    ('dirloop_sse', dict(WHICH_LOOP_TYPE='locopt', **{'MEASURE[Green Function]': True})),
    ('dirloop_sse', dict(LATTICE='triangular lattice', L=3)),
    ('dirloop_sse', dict(MODEL='boson Hubbard', Nmax=2, U=1., t=.3, mu=.5,
                         **{'MEASURE[Local Compressibility]': True,
                            'MEASURE[Site Compressibility]': True,
                            'MEASURE_CORRELATIONS[nn]': 'n:n',
                            'MEASURE_STRUCTURE_FACTOR[nq]': 'n:n'})),
]


@pytest.mark.parametrize('app,parameters', CASES)
@pytest.mark.parametrize('rng', ['mt19937', 'lagged_fibonacci607'])
def test_exact_continuation(executables, tmp_path, app, parameters, rng):
    exe = executables[app]
    run(exe, tmp_path, 'full', parameters, rng=rng)
    for phase, budget in [('thermal', 31), ('production', 247)]:
        run(exe, tmp_path, phase, parameters, rng=rng, budget=budget)
        run(exe, tmp_path, 'resumed', parameters, rng=rng, checkpoint=phase + '.checkpoint.h5')
        compare(tmp_path / 'full.h5', tmp_path / 'resumed.h5')
        compare(tmp_path / 'full.checkpoint.h5', tmp_path / 'resumed.checkpoint.h5')
    run(exe, tmp_path, 'longer', dict(parameters, SWEEPS=1400), rng=rng)
    run(exe, tmp_path, 'extended', dict(parameters, SWEEPS=1400), rng=rng,
        checkpoint='full.checkpoint.h5')
    compare(tmp_path / 'longer.h5', tmp_path / 'extended.h5')


@pytest.mark.parametrize('app', ['worm', 'dirloop_sse'])
@pytest.mark.parametrize('fault', ['chain', 'hamiltonian', 'measurements', 'physical'])
def test_failed_load_preserves_outputs(executables, tmp_path, app, fault):
    exe = executables[app]
    run(exe, tmp_path, 'partial', {}, budget=247)
    with h5py.File(tmp_path / 'partial.checkpoint.h5', 'a') as ar:
        root = ar['simulation/realizations/0/clones/1']
        if fault == 'measurements':
            del root['measurements/Energy']
        else:
            path = fault if fault != 'physical' else ('spins' if app == 'dirloop_sse' else 'initial_state')
            root['checkpoint/' + path][...] = 255
    config = write_run_file(tmp_path / 'bad.toml', parameters=dict(
        LATTICE='chain lattice', MODEL='spin', L=4, local_S=.5, J=1., T=1.,
        THERMALIZATION=100, SWEEPS=1000, SKIP=3, **({'INITIAL_CUTOFF':64} if app=='dirloop_sse' else {})),
        execution=dict(seed=137, bins=8, chains=2), input=dict(checkpoint='partial.checkpoint.h5'),
        output=dict(results='keep.h5', checkpoint='keep.checkpoint.h5'))
    for name in ('keep.h5', 'keep.checkpoint.h5'):
        (tmp_path / name).write_bytes(b'original')
    before = {p.name:p.read_bytes() for p in tmp_path.iterdir()}
    assert subprocess.run([exe, str(config)], capture_output=True, timeout=30).returncode != 0
    assert before == {p.name:p.read_bytes() for p in tmp_path.iterdir()}


def test_chain_density_square_uses_measured_density(executables, tmp_path):
    output = run(executables['worm'], tmp_path, 'chains', dict(
        LATTICE='square lattice', L=2, MODEL='boson Hubbard', Nmax=2,
        U=1., t=.3, mu=.5, CHAIN_KAPPA=True))
    with h5py.File(output) as ar:
        mean = ar['simulation/results/Chain density/mean/value'][()]
        square = ar['simulation/results/Chain density^2/mean/value'][()]
        assert np.all(square >= mean**2 - 1e-12)
        assert np.all(square > 0)


@pytest.mark.parametrize('app,model', [(app, model) for app in ['worm', 'dirloop_sse']
                                     for model in ['spin', 'boson', 'signed']] + [('worm', 'canonical')])
def test_matches_exact_diagonalization(executables, tmp_path, app, model):
    from itertools import product
    from pyalps import alea, hdf5
    sites = 3 if model == 'signed' else 4
    levels = 3 if model in ('boson', 'canonical') else 2
    basis = list(product(range(levels), repeat=sites))
    indices = {state:i for i,state in enumerate(basis)}
    h = np.zeros((len(basis), len(basis)))
    number = np.array([sum(s) for s in basis])
    for row, state in enumerate(basis):
        if model in ('boson', 'canonical'):
            h[row,row] += sum(.5*n*(n-1)-.5*n for n in state)
        for left in range(sites):
            right = (left+1)%sites
            a,b = state[left],state[right]
            if model not in ('boson', 'canonical'):
                h[row,row] += (a-.5)*(b-.5)
                if a != b:
                    target = list(state); target[left],target[right] = b,a
                    h[row,indices[tuple(target)]] += .5
            else:
                for source,destination in ((left,right),(right,left)):
                    if state[source] and state[destination]+1<levels:
                        target=list(state); target[source]-=1; target[destination]+=1
                        h[row,indices[tuple(target)]] -= .3*np.sqrt(state[source]*(state[destination]+1))
    if model == 'canonical':
        selected = number == 4
        h = h[np.ix_(selected, selected)]
        number = number[selected]
    energies, vectors = np.linalg.eigh(h)
    weights = np.exp(-energies); weights /= weights.sum()
    expected = {'Energy Density': np.dot(energies,weights)/sites}
    if model in ('boson', 'canonical'):
        n = np.dot(np.sum(vectors**2*number[:,None],axis=0),weights)
        n2 = np.dot(np.sum(vectors**2*number[:,None]**2,axis=0),weights)
        expected.update(Density=n/sites, Compressibility=(n2-n*n)/sites)
    p=dict(MODEL='spin', local_S=.5, J=1., T=1., THERMALIZATION=2000, SWEEPS=50000)
    p.update(GRAPH='triangle') if model=='signed' else p.update(LATTICE='chain lattice',L=4)
    if model in ('boson', 'canonical'):p.update(MODEL='boson Hubbard',Nmax=2,U=1.,t=.3,mu=.5,NONLOCAL=False)
    if model == 'canonical': p.update(NUMBER_OF_PARTICLES=4, CORRECTION=.1, THERMALIZATION=100000)
    if app=='dirloop_sse':p['INITIAL_CUTOFF']=64
    config=write_run_file(tmp_path/'physics.toml',parameters=p,execution=dict(seed=137,bins=64,chains=3),
                          output=dict(results='physics.h5'))
    completed=subprocess.run([executables[app],str(config)],capture_output=True,text=True,timeout=120)
    assert completed.returncode==0,completed.stdout+completed.stderr
    with hdf5.archive(str(tmp_path/'physics.h5')) as ar:
        for name,value in expected.items():
            result=alea.BatchResult.read(ar,'/simulation/results/'+name)
            assert abs(result.mean[0]-value)<max(.025,5*result.error[0]),(name,result.mean,result.error,value)
    if app=='worm' and model in ('boson', 'canonical'):
        evaluator=str(Path(executables[app]).with_name('worm_evaluate'))
        before=(tmp_path/'physics.h5').read_bytes()
        assert subprocess.run([evaluator,str(tmp_path/'physics.h5')],capture_output=True).returncode==0
        with h5py.File(tmp_path/'physics.h5') as ar:
            assert 'Compressibility' in ar['simulation/results']
            assert 'Centered Density Moments' in ar['simulation/realizations/0/clones/0/results']


@pytest.mark.parametrize('model', ['spin', 'boson', 'signed'])
@pytest.mark.parametrize('skip', [1, 3])
@pytest.mark.parametrize('loop_type', ['minbounce', 'heatbath', 'locopt', 'unweighted'])
def test_green_function_matches_exact_diagonalization(executables, tmp_path, model, skip, loop_type):
    from itertools import product
    from pyalps import alea, hdf5
    sites = 3 if model == 'signed' else 4
    levels = 3 if model == 'boson' else 2
    basis = list(product(range(levels), repeat=sites))
    indices = {state:i for i,state in enumerate(basis)}
    lower = []
    for site in range(sites):
        op = np.zeros((len(basis), len(basis)))
        for col, state in enumerate(basis):
            if state[site]:
                target = list(state); target[site] -= 1
                op[indices[tuple(target)], col] = np.sqrt(state[site])
        lower.append(op)
    number = [op.T @ op for op in lower]
    h = np.zeros_like(lower[0])
    for site in range(sites):
        other = (site + 1) % sites
        hopping = lower[site].T @ lower[other]
        if model == 'boson':
            h += -.3 * (hopping + hopping.T) + .5 * (number[site] @ number[site]) - number[site]
        else:
            z = [n - .5*np.eye(len(basis)) for n in number]
            h += .5 * (hopping + hopping.T) + z[site] @ z[other]
    energy, vectors = np.linalg.eigh(h)
    weights = np.exp(-energy); weights /= weights.sum()
    rho = (vectors * weights) @ vectors.T
    expected = [np.trace(rho @ (.5*(lower[0].T @ op + lower[0] @ op.T))) for op in lower]
    p = dict(MODEL='spin', local_S=.5, J=1., T=1., THERMALIZATION=2000,
             SWEEPS=50000, SKIP=skip, INITIAL_SITE=0, INITIAL_CUTOFF=64,
             **{'MEASURE[Green Function]': True})
    p.update(WHICH_LOOP_TYPE='minbounce' if loop_type == 'unweighted' else loop_type,
             NO_WORMWEIGHT=loop_type == 'unweighted')
    p.update(GRAPH='triangle') if model == 'signed' else p.update(LATTICE='chain lattice', L=4)
    if model == 'boson': p.update(MODEL='boson Hubbard', Nmax=2, U=1., t=.3, mu=.5)
    config = write_run_file(tmp_path/'green.toml', parameters=p,
        execution=dict(seed=137,bins=64,chains=3), output=dict(results='green.h5'))
    completed = subprocess.run([executables['dirloop_sse'], str(config)], capture_output=True, text=True, timeout=120)
    assert completed.returncode == 0, completed.stdout + completed.stderr
    with hdf5.archive(str(tmp_path/'green.h5')) as ar:
        result = alea.BatchResult.read(ar, "/simulation/results/Green's Function")
        assert np.all(np.abs(result.mean-expected) < np.maximum(.02,5*result.error)), (result.mean,result.error,expected)


@pytest.mark.parametrize('target', [0, 4])
def test_canonical_worm_continuation(executables, tmp_path, target):
    p = dict(MODEL='boson Hubbard', Nmax=2, U=1., t=.3, mu=.5,
             NUMBER_OF_PARTICLES=target, CORRECTION=.1, THERMALIZATION=100000, SWEEPS=2000)
    exe = executables['worm']
    run(exe, tmp_path, 'full', p)
    for budget in [31, 147, 100147]:
        run(exe, tmp_path, 'partial', p, budget=budget)
        run(exe, tmp_path, 'resumed', p, checkpoint='partial.checkpoint.h5')
        compare(tmp_path/'full.h5', tmp_path/'resumed.h5')
        compare(tmp_path/'full.checkpoint.h5', tmp_path/'resumed.checkpoint.h5')
    with h5py.File(tmp_path/'full.h5') as ar:
        assert ar['simulation/results/Density/mean/value'][0] == target/4
        assert abs(ar['simulation/results/Compressibility/mean/value'][0]) < 1e-12
        if target == 0:
            assert abs(ar['simulation/results/Energy/mean/value'][0]) < 1e-12


@pytest.mark.parametrize('app,extra,reason', [
    ('dirloop_sse', {'MEASURE[Green Function]':True, 'WORM_ABORT':1.}, 'WORM_ABORT'),
    ('dirloop_sse', {'MEASURE[Green Function]':True, 'RESTRICT_MEASUREMENTS[N]':4}, 'restricted'),
    ('worm', {'NUMBER_OF_PARTICLES':9, 'CORRECTION':.1}, 'allowed range'),
    ('worm', {'NUMBER_OF_PARTICLES':4, 'CORRECTION':0.}, 'positive'),
    ('worm', {'NUMBER_OF_PARTICLES':4, 'CORRECTION':.1, 'ADJUST':'U'}, 'particle-number sector'),
])
def test_invalid_estimator_configuration(executables, tmp_path, app, extra, reason):
    p = dict(LATTICE='chain lattice', L=4, MODEL='boson Hubbard', Nmax=2,
             U=1., t=.3, mu=.5, T=1., THERMALIZATION=100000, SWEEPS=100)
    p.update(extra)
    config = write_run_file(tmp_path/'invalid.toml', parameters=p,
        execution=dict(seed=137), output=dict(results='invalid.h5'))
    done = subprocess.run([executables[app], str(config)], capture_output=True, text=True, timeout=30)
    assert done.returncode != 0
    assert reason in done.stderr
    assert not (tmp_path/'invalid.h5').exists()


def test_chain_density_on_unequal_fractional_coordinate_groups(executables, tmp_path):
    vertices = ''.join(f'<VERTEX id="{i+1}"><COORDINATE>{i} {y}</COORDINATE></VERTEX>'
                       for i,y in enumerate([-.5,-.5,1.25,1.25,1.25]))
    edges = ''.join(f'<EDGE source="{i+1}" target="{(i+1)%5+1}"/>' for i in range(5))
    (tmp_path/'graph.xml').write_text('<LATTICES><GRAPH name="chains" vertices="5" dimension="2">'
                                    +vertices+edges+'</GRAPH></LATTICES>')
    config = write_run_file(tmp_path/'chains.toml', parameters=dict(
        GRAPH='chains', MODEL='boson Hubbard', Nmax=2, U=1., t=.3, mu=.5, T=1.,
        THERMALIZATION=100, SWEEPS=1000, CHAIN_KAPPA=True),
        input=dict(lattice_library='graph.xml'), execution=dict(seed=137,bins=8),
        output=dict(results='chains.h5'))
    done = subprocess.run([executables['worm'],str(config)],capture_output=True,text=True,timeout=30)
    assert done.returncode == 0, done.stdout+done.stderr
    with h5py.File(tmp_path/'chains.h5') as ar:
        results = ar['simulation/results']
        mean = results['Chain density/mean/value'][()]
        square = results['Chain density^2/mean/value'][()]
        assert mean.shape == (2,)
        assert np.all(square >= mean**2 - 1e-12)
        np.testing.assert_allclose(np.dot([2,3],mean)/5, results['Density/mean/value'][0])
        assert 'Stiffness' not in results  # Explicit graph has no periodic bond vectors.
