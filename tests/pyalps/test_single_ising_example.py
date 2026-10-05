"""Heat-bath Ising and warm-start scans preserve physics and exact continuation."""
import itertools
import os
from pathlib import Path
import subprocess

import h5py
import numpy as np
import pyalps
import pytest
from pyalps import alea, hdf5
from pyalps.run_io import write_run_file
from test_native_mpi import compare, invoke, launcher


@pytest.fixture
def executable():
    root=os.environ.get('ALPS_EXAMPLES_BUILD_DIR')
    if not root:
        pytest.skip('Set ALPS_EXAMPLES_BUILD_DIR to the built C++ examples')
    return str(Path(root).resolve()/'parapack/single/ising_single')


def run(executable, directory, name, *, scan=False, budget=0, checkpoint=None,
        rng='mt19937', sweeps=117, coupling=1., lattice='chain lattice', length=4,
        processes=1, launcher=None, threads=None):
    parameters=dict(LATTICE=lattice,L=length,J=coupling,SWEEPS=sweeps,THERMALIZATION=11)
    if scan:
        parameters.update(ALGORITHM='ising; temperature scan',NUM_TEMPERATURES=3,
                          INITIAL_TEMPERATURE=2.,DIFF_TEMPERATURE=-.5,INITIAL_THERMALIZATION=19)
    else:
        parameters['T']=2.
    path=directory/(name+'.toml')
    write_run_file(path, parameters=parameters,
                   input={'checkpoint':str(checkpoint)} if checkpoint else {},
                   output=dict(results=name+'.h5',checkpoint=name+'.checkpoint.h5'),
                   execution=dict(seed=823,chains=2,bins=8,rng=rng,max_sweeps=budget))
    if threads:
        environment=dict(os.environ,OMP_NUM_THREADS=str(threads))
        process=subprocess.run([executable,str(path)],capture_output=True,text=True,env=environment)
        assert process.returncode==0,process.stdout+process.stderr
    else:
        invoke(launcher,executable,path,processes=processes)
    return directory/(name+'.h5'),directory/(name+'.checkpoint.h5')


@pytest.mark.parametrize('rng',['mt19937','lagged_fibonacci607'])
@pytest.mark.parametrize('scan,budget',[(False,5),(False,38),(True,9),(True,71),(True,142),(True,301)])
def test_exact_continuation(executable,tmp_path,rng,scan,budget):
    full,full_state=run(executable,tmp_path,'full',rng=rng,scan=scan)
    _,state=run(executable,tmp_path,'part',rng=rng,scan=scan,budget=budget)
    resumed,resumed_state=run(executable,tmp_path,'resumed',rng=rng,scan=scan,checkpoint=state)
    compare(full,resumed);compare(full_state,resumed_state)
    # Different thread counts must consume the same RNG stream and updates.
    threaded,threaded_state=run(executable,tmp_path,'threaded',rng=rng,scan=scan,threads=2)
    compare(full,threaded);compare(full_state,threaded_state)


@pytest.mark.parametrize('fault',['spin','count','topology','rng'])
def test_extension_and_corrupt_checkpoint(executable,tmp_path,fault):
    full,full_state=run(executable,tmp_path,'full')
    _,state=run(executable,tmp_path,'short',sweeps=23)
    extended,extended_state=run(executable,tmp_path,'extended',checkpoint=state)
    compare(full,extended);compare(full_state,extended_state)
    with h5py.File(state,'r+') as archive:
        state=archive['simulation/realizations/0/clones/1/checkpoint']
        if fault=='spin': state['spins'][0]=0
        elif fault=='count': state['sweeps'][()]+=1
        elif fault=='topology': state['topology'][1]=999
        else: del state['engine']
    before=extended.read_bytes(),extended_state.read_bytes()
    process=subprocess.run([executable,str(tmp_path/'extended.toml')],capture_output=True)
    assert process.returncode!=0
    assert before==(extended.read_bytes(),extended_state.read_bytes())


@pytest.mark.parametrize('coupling',[1.,-1.])
def test_scan_exact_physics_and_analysis(executable,tmp_path,coupling):
    filename,_=run(executable,tmp_path,'physics',scan=True,sweeps=50000,coupling=coupling)
    datasets=pyalps.loadMeasurements([str(filename)])
    assert len(datasets)==3
    spins=np.array(list(itertools.product((-1.,1.),repeat=4)))
    energy=-coupling*np.sum(spins*np.roll(spins,1,axis=1),axis=1)
    mag=spins.sum(axis=1)
    values=np.array([np.full(16,4.),energy,energy**2,mag,mag**2,mag**4]).T
    for i,(temperature,data) in enumerate(zip((2.,1.5,1.),datasets)):
        expected=np.average(values,axis=0,weights=np.exp(-energy/temperature))
        names=('Number of Sites','Energy','Energy^2','Magnetization','Magnetization^2','Magnetization^4')
        results={d.props['observable']:d.native_result for d in data}
        assert all(d.props['T']==temperature for d in data)
        for j,name in enumerate(names):
            assert results[name].count==100000
            np.testing.assert_allclose(results[name].mean,[expected[j]],atol=max(.015*max(1,abs(expected[j])),6*results[name].error[0]))
        path=f'/simulation/replicas/{i}'
        with hdf5.archive(str(filename)) as archive:
            joint=alea.read_result(archive,path+'/joint')
        keep=joint.batch_counts>0
        weights,sums=joint.batch_counts[keep].astype(float),joint.batch_sums[keep]
        count=weights.sum();total=sums.sum(axis=0);mean=total/count
        leave=(total-sums)/(count-weights[:,None])
        pseudo=count*mean[4]**2/mean[5]-(count-weights)*leave[:,4]**2/leave[:,5]
        np.testing.assert_allclose(results['Binder Ratio of Magnetization'].batch_sums[keep,0],pseudo,rtol=1e-10)
        heat_pseudo=(count*(mean[2]-mean[1]**2)
                     -(count-weights)*(leave[:,2]-leave[:,1]**2))/(temperature**2*4.)
        np.testing.assert_allclose(results['Specific Heat'].batch_sums[keep,0],heat_pseudo,rtol=1e-10,atol=1e-8)
        exact_heat=(expected[2]-expected[1]**2)/(temperature**2*4.)
        np.testing.assert_allclose(results['Specific Heat'].mean,[exact_heat],atol=max(.015,6*results['Specific Heat'].error[0]))
    diagnostics=pyalps.loadBinningAnalysis([str(filename)])
    assert len(diagnostics)==3 and all(len(group)==1 and group[0].native_result.count==50000 for group in diagnostics)


def test_mpi_scan_checkpoint_repartition(executable,launcher,tmp_path):
    full,full_state=run(executable,tmp_path,'full',scan=True)
    _,state=run(executable,tmp_path,'part',scan=True,budget=142,processes=2,launcher=launcher)
    resumed,resumed_state=run(executable,tmp_path,'resumed',scan=True,checkpoint=state)
    compare(full,resumed);compare(full_state,resumed_state)


def test_thread_count_independent_large_lattice(executable,tmp_path):
    one,one_state=run(executable,tmp_path,'one',length=512,sweeps=37,scan=True,threads=1)
    two,two_state=run(executable,tmp_path,'two',length=512,sweeps=37,scan=True,threads=2)
    four,four_state=run(executable,tmp_path,'four',length=512,sweeps=37,scan=True,threads=4)
    compare(one,two);compare(one,four)
    compare(one_state,two_state);compare(one_state,four_state)


def test_uncoupled_extreme_beta_and_validation(executable,tmp_path):
    path=tmp_path/'uncoupled.toml'
    write_run_file(path,parameters=dict(LATTICE='chain lattice',L=4,J=0.,T=1.e-300,
                                        SWEEPS=20000,THERMALIZATION=0),
                   execution=dict(seed=823,bins=32),output=dict(results='uncoupled.h5'))
    process=subprocess.run([executable,str(path)],capture_output=True)
    assert process.returncode==0,process.stderr
    data={d.props['observable']:d.native_result for d in pyalps.loadMeasurements([str(tmp_path/'uncoupled.h5')])[0]}
    np.testing.assert_array_equal(data['Energy'].mean,[0.])
    np.testing.assert_array_equal(data['Specific Heat'].mean,[0.])
    np.testing.assert_allclose(data['Magnetization^2'].mean,[4.],atol=.12)
    np.testing.assert_allclose(data['Magnetization^4'].mean,[40.],atol=1.5)
    original=(tmp_path/'uncoupled.h5').read_bytes()
    write_run_file(path,parameters=dict(LATTICE='chain lattice',L=4,J=0.,T=0.,
                                        SWEEPS=20000,THERMALIZATION=0),
                   execution=dict(seed=823,bins=32),output=dict(results='uncoupled.h5'),overwrite=True)
    assert subprocess.run([executable,'--validate',str(path)],capture_output=True).returncode!=0
    assert (tmp_path/'uncoupled.h5').read_bytes()==original


def test_scan_first_stage_matches_fixed_temperature(executable,tmp_path):
    scanned,scan_state=run(executable,tmp_path,'scan',scan=True,budget=136)
    path=tmp_path/'fixed.toml'
    write_run_file(path,parameters=dict(LATTICE='chain lattice',L=4,J=1.,T=2.,
                                        SWEEPS=117,THERMALIZATION=19),
                   execution=dict(seed=823,chains=2,bins=8),
                   output=dict(results='fixed.h5',checkpoint='fixed.checkpoint.h5'))
    invoke(None,executable,path)
    with hdf5.archive(str(scanned)) as archive:
        scan_joint=alea.read_result(archive,'/simulation/replicas/0/joint')
    with hdf5.archive(str(tmp_path/'fixed.h5')) as archive:
        fixed_joint=alea.read_result(archive,'/simulation/joint')
    np.testing.assert_array_equal(scan_joint.batch_sums,fixed_joint.batch_sums)
    np.testing.assert_array_equal(scan_joint.batch_counts,fixed_joint.batch_counts)
    with h5py.File(scan_state) as scan,h5py.File(tmp_path/'fixed.checkpoint.h5') as fixed:
        for i in range(2):
            path=f'simulation/realizations/0/clones/{i}/checkpoint/spins'
            np.testing.assert_array_equal(scan[path][()],fixed[path][()])
