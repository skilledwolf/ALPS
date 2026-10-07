"""Energy Wang-Landau: exact DOS, microcanonical reweighting and continuation."""
import itertools
from pathlib import Path
import subprocess

import h5py
import numpy as np
import pyalps
import pytest
from pyalps.run_io import write_run_file
from conftest import tutorials_build
from test_native_mpi import compare, invoke, launcher

BONDS = [(0, 1), (1, 2), (2, 3), (3, 0), (1, 1)]
SPINS = np.array(list(itertools.product((-1, 1), repeat=4)))


@pytest.fixture
def executable():
    return str(tutorials_build()/'00-examples/mc/wanglandau/wanglandau')


def spectrum(coupling=1):
    energies = -coupling*sum(SPINS[:, a]*SPINS[:, b] for a, b in BONDS)
    return energies, SPINS.mean(axis=1)


def weights_file(path, coupling=1, window=None, offset=0.):
    energies, _ = spectrum(coupling)
    low, high = window or (int(energies.min()), int(energies.max()))
    grid = np.arange(low, high+1)
    counts = (grid[:, None] == energies).sum(axis=1)
    with h5py.File(path, 'w') as ar:
        g = ar.create_group('weights')
        g['range'] = np.array([low, high], dtype=np.int64)
        g['sites'] = np.uint64(4)
        g['coupling'] = np.int64(coupling)
        g['bonds'] = np.array(BONDS, dtype=np.uint64).ravel()
        g['complete'] = True
        g['visited'] = counts > 0
        g['logg'] = np.log(np.maximum(counts, 1)) + offset
    return path


def run(directory, executable, name, *, mode='learn', coupling=1, parameters=None,
        input=None, execution=None, checkpoint=True, processes=1, launcher=None):
    energy, _ = spectrum(coupling)
    p = dict(MODE=mode, GRAPH='test', COUPLING=coupling,
             ENERGY_RANGE=[int(energy.min()), int(energy.max())], SWEEPS=397,
             THERMALIZATION=31, CHECK_INTERVAL=256, FINAL_UPDATE_FACTOR=float(np.exp(.0001)))
    p.update(parameters or {})
    graph = directory/'graph.xml'
    graph.write_text('<LATTICES><GRAPH name="test" vertices="4" dimension="1">'+''.join(
        f'<EDGE source="{a+1}" target="{b+1}"/>' for a, b in BONDS)+'</GRAPH></LATTICES>')
    source = dict(lattice_library=str(graph), **(input or {}))
    e = dict(seed=718, bins=16, chains=1)
    e.update(execution or {})
    destination = dict(results=name+'.h5')
    if checkpoint and mode!='reweight': destination['checkpoint'] = name+'.checkpoint.h5'
    path = write_run_file(directory/(name+'.toml'), parameters=p, input=source,
                          execution=e, output=destination, overwrite=True)
    invoke(launcher, executable, path, processes=processes)
    return directory/(name+'.h5'), directory/(name+'.checkpoint.h5')


# The free spins have a single energy, which is both ends of the flat range.
@pytest.mark.parametrize('coupling', [-1, 0])
def test_learning_matches_exact_density(executable, tmp_path, coupling):
    result, state = run(tmp_path, executable, 'learn', coupling=coupling)
    energies, _ = spectrum(coupling)
    with h5py.File(result) as ar:
        g = ar['weights']
        assert g['complete'][()]
        grid = np.arange(energies.min(), energies.max()+1)
        count = (grid[:, None] == energies).sum(axis=1)
        np.testing.assert_array_equal(g['visited'][()], count > 0)
        logg = g['logg'][()][count > 0]
        expected = np.log(count[count > 0])
        np.testing.assert_allclose(logg-logg[0], expected-expected[0], atol=.11)
        clone = ar['simulation/realizations/0/clones/0']
        # Both energy endpoints participate in flatness, including a singleton.
        hist = clone['Final Histogram'][()][count > 0]
        assert hist.min() >= .95*hist.mean()
        assert clone['Overall Histogram'][()].sum() == 4*clone['completed_sweeps'][()]


@pytest.mark.parametrize('mode,rng', [('learn', 'mt19937'), ('measure', 'lagged_fibonacci607')])
def test_exact_restart_across_refinements_and_partial_bins(executable, tmp_path, rng, mode):
    source = {} if mode=='learn' else dict(weights=[str(weights_file(tmp_path/'weights.h5'))])
    options = dict(mode=mode, input=source, execution=dict(rng=rng))
    full, full_state = run(tmp_path, executable, 'full', **options)
    for budget in ((7, 319, 617) if mode == 'learn' else (7, 319)):
        _, state = run(tmp_path, executable, f'part{budget}', **dict(options, execution=dict(rng=rng, max_sweeps=budget)))
        result, restored = run(tmp_path, executable, f'resume{budget}', **dict(options, input=dict(source, checkpoint=str(state))))
        compare(full, result); compare(full_state, restored)
    if mode=='measure':
        _, reference = run(tmp_path, executable, 'longer', parameters=dict(SWEEPS=677), **options)
        _, extended = run(tmp_path, executable, 'extended', parameters=dict(SWEEPS=677), **dict(options, input=dict(source, checkpoint=str(full_state))))
        compare(reference, extended)


def test_empty_measurement_interval_cannot_refine(executable, tmp_path):
    previous_stage,previous_visits=0,0
    empty_intervals=0
    # Two refinements, each followed by intervals that never reach the measured energy.
    for budget in range(1,8):
        _,state=run(tmp_path,executable,'trace',
            parameters=dict(CHECK_INTERVAL=1,ENERGY_MEASURE_RANGE=[-5,-5],
                            VISIT_PENALTY=1.,INITIAL_UPDATE_FACTOR=1.01),
            execution=dict(max_sweeps=budget))
        with h5py.File(state) as ar:
            clone=ar['simulation/realizations/0/clones/0']
            stage,visits=int(clone['stage'][()]),int(clone['overall'][0])
        if visits==previous_visits and previous_visits:
            empty_intervals+=1
            assert stage==previous_stage
        previous_stage,previous_visits=stage,visits
    assert empty_intervals


@pytest.mark.parametrize('coupling', [1])
def test_stitched_windows_and_canonical_physics(executable, tmp_path, coupling):
    energy, mag = spectrum(coupling)
    low, high = int(energy.min()), int(energy.max())
    middle = sorted(set(energy))[len(set(energy))//2]
    left = weights_file(tmp_path/'left.h5', coupling, [low, int(middle)], offset=-500.)
    right = weights_file(tmp_path/'right.h5', coupling, [int(middle), high], offset=700.)
    result, _ = run(tmp_path, executable, 'measure', mode='measure', coupling=coupling,
                    parameters=dict(SWEEPS=100000), input=dict(weights=[str(right), str(left)]),
                    execution=dict(chains=3))
    temperatures = [.6, 1.2, 8.]
    analysis, _ = run(tmp_path, executable, 'analysis', mode='reweight', coupling=coupling,
        parameters=dict(TEMPERATURE_SET=temperatures, REFERENCE_BIN=low,
                        REFERENCE_LOGG=float(np.log(np.count_nonzero(energy==low)))),
        input=dict(measurements=[str(result)]))
    groups = pyalps.loadMeasurements([str(analysis)])
    assert len(groups)==len(temperatures)
    for temperature, group in zip(temperatures, groups):
        beta=1/temperature
        boltzmann=np.exp(-beta*energy)
        average=lambda x: np.average(x, weights=boltzmann)
        e,m2,m4=average(energy),average(mag**2),average(mag**4)
        targets={'Energy': e, 'Energy Density': e/4, 'Energy^2': average(energy**2),
                 'Specific Heat': beta**2*(average(energy**2)-e**2)/4,
                 'Magnetization': average(mag), 'Magnetization^2': m2, 'Magnetization^4': m4,
                 'Binder Ratio of Magnetization': m2*m2/m4,
                 'Free-Energy': -np.log(boltzmann.sum())/beta,
                 'Entropy': np.log(boltzmann.sum())+beta*e}
        results={d.props['observable']: d.native_result for d in group}
        for name, target in targets.items():
            actual=results[name]
            assert actual.count==300000
            np.testing.assert_allclose(actual.mean, [target], atol=max(.008,7*actual.error[0]), err_msg=name)
        assert group[0].props['T']==temperature


def test_distinct_walk_measure_windows_and_unavailable_normalization(executable, tmp_path):
    weights=weights_file(tmp_path/'window.h5', window=[-5, -1])
    result, _ = run(tmp_path, executable, 'measure', mode='measure',
        parameters=dict(ENERGY_MEASURE_RANGE=[-5, -1], SWEEPS=1000), input=dict(weights=[str(weights)]))
    with h5py.File(result) as ar:
        clone=ar['simulation/realizations/0/clones/0']
        assert clone['Overall Histogram'][-1]>0
        assert clone['Measurement Histogram'][-1]==0
    analysis, _ = run(tmp_path, executable, 'analysis', mode='reweight',
        parameters=dict(ENERGY_MEASURE_RANGE=[-5, -1], TEMPERATURE_SET=[1.]),
        input=dict(measurements=[str(result)]))
    with h5py.File(analysis) as ar:
        assert 'Energy' in ar['simulation/replicas/0/results']
        assert 'Free-Energy' in ar['simulation/replicas/0/unavailable']
        assert 'Entropy' in ar['simulation/replicas/0/unavailable']


def test_corrupt_checkpoint_does_not_replace_outputs(executable, tmp_path):
    weights=weights_file(tmp_path/'weights.h5')
    options=dict(mode='measure', input=dict(weights=[str(weights)]))
    _, state=run(tmp_path, executable, 'partial', execution=dict(max_sweeps=67), **options)
    result, resumed=run(tmp_path, executable, 'resume', **dict(options, input=dict(options['input'], checkpoint=str(state))))
    before=result.read_bytes(),resumed.read_bytes()
    original = state.read_bytes()
    for fault in ['spin', 'energy', 'topology', 'logg', 'count', 'histogram', 'rng']:
        state.write_bytes(original)
        with h5py.File(state, 'r+') as ar:
            clone=ar['simulation/realizations/0/clones/0']
            if fault=='spin': clone['spins'][0]=0
            elif fault=='energy': clone['energy'][()]+=1
            elif fault=='topology': clone['weights/bonds'][0]=999
            elif fault=='logg': clone['weights/logg'][0]+=1
            elif fault=='count': clone['steps'][()]+=1
            elif fault=='histogram': clone['overall'][0]=2**63
            else: del clone['rng']
        rejected=subprocess.run([executable,str(tmp_path/'resume.toml')],capture_output=True,timeout=30)
        assert rejected.returncode!=0
        assert before==(result.read_bytes(),resumed.read_bytes())


@pytest.mark.parametrize('mode', ['learn', 'measure'])
def test_mpi_independent_windows_and_cross_rank_restart(executable, tmp_path, launcher, mode):
    source={} if mode=='learn' else dict(weights=[str(weights_file(tmp_path/'weights.h5'))])
    options=dict(mode=mode, input=source, execution=dict(chains=3))
    full, full_state=run(tmp_path, executable, 'full', **options)
    _, state=run(tmp_path, executable, 'partial', **dict(options, execution=dict(chains=3,max_sweeps=127)), processes=2,launcher=launcher)
    result, restored=run(tmp_path, executable, 'resume', **dict(options,input=dict(source,checkpoint=str(state))), processes=3,launcher=launcher)
    compare(full,result);compare(full_state,restored)


def test_learning_penalty_and_native_pipeline(executable, tmp_path):
    partial, state=run(tmp_path, executable, 'penalty', execution=dict(max_sweeps=31),
                      parameters=dict(ENERGY_MEASURE_RANGE=[-5, -1], VISIT_PENALTY=5.))
    with h5py.File(state) as ar:
        clone=ar['simulation/realizations/0/clones/0']
        counts=clone['overall'][()]
        penalty=np.where(np.arange(-5,4)<=-1,1.,5.)
        np.testing.assert_allclose(clone['weights/logg'][()],counts*penalty)
        assert counts[-1]>0
    learned, _=run(tmp_path, executable, 'learned')
    observations=[]
    for i in range(2):
        result, _=run(tmp_path, executable, f'measured{i}', mode='measure',
                      execution=dict(seed=124+i), input=dict(weights=[str(learned)]))
        observations.append(str(result))
    analysis, _=run(tmp_path, executable, 'analysis', mode='reweight',
                    parameters=dict(TEMPERATURE_SET=[2.], REFERENCE_BIN=-5, REFERENCE_LOGG=float(np.log(2.))),
                    input=dict(measurements=observations))
    values={d.props['observable']:d.native_result for d in pyalps.loadMeasurements([str(analysis)])[0]}
    energy, _=spectrum()
    assert values['Energy'].count==2*397
    np.testing.assert_allclose(values['Energy'].mean,[np.average(energy,weights=np.exp(-energy/2))],
                               atol=max(.1,7*values['Energy'].error[0]))


def test_validation_is_transactional(executable, tmp_path):
    weights=weights_file(tmp_path/'weights.h5')
    result, _=run(tmp_path, executable, 'base', mode='measure', input=dict(weights=[str(weights)]))
    import tomllib
    weights_data, result_data = weights.read_bytes(), result.read_bytes()
    for fault in ('missing-file', 'empty-files', 'duplicate-files', 'output-input',
                  'incomplete-weights', 'incompatible-model', 'disjoint-windows', 'inconsistent-coverage',
                  'nonfinite-weight', 'invalid-range', 'wrong-measure-range', 'missing-reference',
                  'incomplete-measurements', 'different-measurement-weights', 'invalid-temperature', 'invalid-factor'):
        weights.write_bytes(weights_data)
        result.write_bytes(result_data)
        document=tomllib.loads((tmp_path/'base.toml').read_text())
        p,source,output=document['parameters'],document['input'],document['output']
        output.update(results='keep.h5',checkpoint='keep.checkpoint.h5')
        if fault=='missing-file': source['weights'].append(str(tmp_path/'absent.h5'))
        elif fault=='empty-files': source['weights']=[]
        elif fault=='duplicate-files': source['weights'].append(str(weights))
        elif fault=='output-input': output['results']=str(weights)
        elif fault in ('incomplete-weights','incompatible-model','nonfinite-weight'):
            with h5py.File(weights,'r+') as ar:
                if fault=='incomplete-weights': ar['weights/complete'][()]=False
                elif fault=='incompatible-model': ar['weights/coupling'][()]=-1
                else: ar['weights/logg'][0]=np.nan
        elif fault in ('disjoint-windows','inconsistent-coverage'):
            left=weights_file(tmp_path/'left.h5',window=[-5,-1])
            right=weights_file(tmp_path/'right.h5',window=[0,3] if fault=='disjoint-windows' else [-1,3])
            if fault=='inconsistent-coverage':
                with h5py.File(right,'r+') as ar: ar['weights/visited'][0]=False
            source['weights']=[str(left),str(right)]
        elif fault=='invalid-range': p['ENERGY_RANGE']=[3,-5]
        elif fault=='wrong-measure-range': p['ENERGY_MEASURE_RANGE']=[-6,3]
        elif fault=='invalid-factor':
            p.update(MODE='learn',FINAL_UPDATE_FACTOR=1.)
            del source['weights']
        else:
            p.update(MODE='reweight',TEMPERATURE_SET=[2.])
            del source['weights'];source['measurements']=[str(result)]
            del output['checkpoint']
            if fault=='missing-reference': p['REFERENCE_BIN']=-5
            elif fault=='invalid-temperature': p['TEMPERATURE_SET']=[0.]
            elif fault=='incomplete-measurements':
                with h5py.File(result,'r+') as ar: ar['microcanonical/complete'][()]=False
            else:
                import shutil
                second=tmp_path/'different.h5';shutil.copyfile(result,second)
                with h5py.File(second,'r+') as ar: ar['weights/logg'][0]+=1.
                source['measurements'].append(str(second))
        for path in ('keep.h5','keep.checkpoint.h5'): (tmp_path/path).write_bytes(b'existing results')
        invalid=write_run_file(tmp_path/'invalid.toml',overwrite=True,**document)
        before={path.name:path.read_bytes() for path in tmp_path.iterdir()}
        rejected=subprocess.run([executable,'--validate',str(invalid)],capture_output=True,timeout=30)
        assert rejected.returncode!=0,(fault,rejected.stdout,rejected.stderr)
        assert before=={path.name:path.read_bytes() for path in tmp_path.iterdir()}


def test_corrupt_learning_state_is_rejected(executable, tmp_path):
    _, state=run(tmp_path,executable,'partial',execution=dict(max_sweeps=17))
    result,resumed=run(tmp_path,executable,'resumed',input=dict(checkpoint=str(state)))
    original = state.read_bytes()
    for fault in ['done', 'stage', 'histogram', 'nonfinite']:
        state.write_bytes(original)
        with h5py.File(state,'r+') as ar:
            clone=ar['simulation/realizations/0/clones/0']
            if fault=='done': clone['done'][()]=True
            elif fault=='stage': clone['stage'][()]=1
            elif fault=='histogram': clone['histogram'][0]=2**63
            else: clone['weights/logg'][0]=np.nan
        before=result.read_bytes(),resumed.read_bytes()
        rejected=subprocess.run([executable,str(tmp_path/'resumed.toml')],capture_output=True,timeout=30)
        assert rejected.returncode!=0
        assert before==(result.read_bytes(),resumed.read_bytes())
