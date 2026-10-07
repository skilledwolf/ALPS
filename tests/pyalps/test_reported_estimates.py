"""Published estimates retain only the statistics that were actually recorded."""
import copy
from pathlib import Path
import subprocess
import sys

import h5py
import numpy as np
import pytest
import pyalps
from pyalps import alea, hdf5
from pyalps.dataset import DataSet
from test_analysis_commands import command

ROOT = Path(__file__).resolve().parents[2]


@pytest.mark.parametrize('complex_mean', [False, True])
@pytest.mark.parametrize('optional', [False, True])
def test_reported_estimate_roundtrip_and_measurement_io(tmp_path, complex_mean, optional):
    mean = np.array([2., 4.]) * (1+2j if complex_mean else 1)
    fields = dict(error=np.array([.3, .7]), variance=np.array([8., 9.]),
                  tau=np.array([1., 2.]), converged_errors=np.array([0, 2])) if optional else {}
    result = alea.ReportedEstimate(count=71, mean=mean, **fields)
    mean[:] = 9  # the record owns its arrays
    assert result.mean[0] != 9
    assert not hasattr(result, 'count2') and not hasattr(result, 'batch_counts')
    assert not hasattr(result, 'covariance')
    filename = tmp_path/'reported.h5'
    with hdf5.archive(filename, 'w') as ar:
        ar['/estimate'] = result
    with hdf5.archive(filename) as ar:
        loaded = alea.read_result(ar, '/estimate')
    assert type(loaded) is alea.ReportedEstimate and loaded.count == 71
    for name in result._fields:
        assert hasattr(loaded, name) == hasattr(result, name)
        if hasattr(result, name):
            np.testing.assert_array_equal(getattr(loaded, name), getattr(result, name))
    data = DataSet.from_result(result, props={'observable': 'Reported / E'})
    pyalps.saveMeasurements([data], filename)
    reloaded = pyalps.loadMeasurements([str(filename)])[0][0]
    np.testing.assert_array_equal(reloaded.native_result.mean, result.mean)
    assert hasattr(reloaded.native_result, 'error') == optional
    with pytest.raises(ValueError, match='reported estimates'):
        pyalps.mergeDataSets([data, data])


def test_reported_save_removes_absent_owned_fields_and_preserves_metadata(tmp_path):
    rich = alea.ReportedEstimate(count=17, mean=[2.], error=[.5], variance=[3.], tau=[1.])
    minimal = alea.ReportedEstimate(count=5, mean=[7.])
    with hdf5.archive(tmp_path/'replace.h5', 'w') as ar:
        rich.save(ar, '/estimate')
        ar['/estimate/variance/note'] = 'kept'
        ar['/estimate/legacy/count'] = np.uint64(17)
        minimal.save(ar, '/estimate')
        assert not ar.is_data('/estimate/variance/value')
        assert not ar.is_data('/estimate/mean/error')
        assert ar['/estimate/variance/note'] == 'kept'
        assert ar['/estimate/legacy/count'] == 17
        invalid = copy.deepcopy(minimal)
        invalid.error = np.array([1., 2.])
        with pytest.raises(ValueError):
            invalid.save(ar, '/estimate')
        assert ar['/estimate/count'] == 5
        np.testing.assert_array_equal(ar['/estimate/mean/value'], [7.])
        alea.MeanAccumulator().result().save(ar, '/native')
        with pytest.raises(ValueError, match='native estimator'):
            minimal.save(ar, '/native')


@pytest.mark.parametrize('fields', [dict(count=-1, mean=[0.]), dict(count=True, mean=[0.]),
    dict(count=2**64, mean=[0.]), dict(count=1, mean=0.), dict(count=1, mean=[]),
    dict(count=1, mean=[1.], error=[1., 2.]), dict(count=1, mean=[1.], count2=1),
    dict(count=1, mean=[1.], error=[1j])])
def test_invalid_reported_records(fields):
    with pytest.raises(ValueError):
        alea.ReportedEstimate(**fields)


@pytest.mark.parametrize('vector,history', [(False, True), (True, False)])
def test_offline_summary_conversion_and_analysis_commands(tmp_path, vector, history):
    source, target = tmp_path/'legacy.h5', tmp_path/'converted.h5'
    mean = np.array([2., 5.]) if vector else np.array(2.)
    with h5py.File(source, 'w') as ar:
        ar.create_group('parameters')
        group = ar.create_group('simulation/results/E')
        group['count'] = np.uint64(81)
        group['mean/value'], group['mean/error'] = mean, mean/10
        group['variance/value'], group['tau/value'] = mean*3, mean/2
        group.attrs['cannotrebin'] = True
        if history:
            group['timeseries/data'] = np.stack([mean, mean*2])
            group['timeseries/data'].attrs.update(binsize=40, binningtype='linear')
    before = source.read_bytes()
    conversion = [sys.executable, str(ROOT/'src/tools/hdf5/convert.py'), str(source), str(target),
                  '--alea-summary', '/simulation/results/E']
    completed = subprocess.run(conversion, text=True, capture_output=True)
    assert completed.returncode == 0, completed.stderr
    assert source.read_bytes() == before
    data = pyalps.loadMeasurements([str(target)])[0][0]
    result = data.native_result
    assert isinstance(result, alea.ReportedEstimate) and result.count == 81
    np.testing.assert_array_equal(result.mean, np.atleast_1d(mean))
    np.testing.assert_array_equal(result.error, np.atleast_1d(mean/10))
    np.testing.assert_array_equal(result.variance, np.atleast_1d(mean*3))
    np.testing.assert_array_equal(result.tau, np.atleast_1d(mean/2))
    with h5py.File(target) as ar:
        group = ar['simulation/results/E']
        assert 'kind' not in group.attrs
        assert 'batch' not in group and 'count2' not in group
        assert ('legacy/timeseries/data' in group) == history
        np.testing.assert_array_equal(group['legacy/mean/value'], mean)
    for tool in ('mean', 'variance'):
        completed = command(tool, '-v', '-w', target)
        assert completed.returncode == 0, completed.stderr


def test_missing_reported_variance_is_not_invented_by_cli(tmp_path):
    path = tmp_path/'mean-only.h5'
    with hdf5.archive(path, 'w') as ar:
        alea.ReportedEstimate(count=12, mean=[2.], error=[.7]).save(ar, '/simulation/results/E')
    before = path.read_bytes()
    completed = command('variance', '-w', path)
    assert completed.returncode != 0 and 'variance cannot be recovered' in completed.stderr
    assert path.read_bytes() == before
