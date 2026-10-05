"""Native result discovery preserves encoded names and component vectors."""
import numpy as np
import pytest

import pyalps
from pyalps import alea, hdf5

RESULTS = '/simulation/results'
VALUES = {
    'Energy': ([-1.2345], [.0067]),
    'Magnetization^2': ([.3721], [.0014]),
    'Energy Density / Site': ([-.6172], [.0033]),
    'Correlations': ([1., .5, .25], [.01, .02, .03]),
    "Green's Function": ([.9, -.4], [.05, .06]),
}


@pytest.fixture
def results_archive(tmp_path):
    path = str(tmp_path/'loadobs.h5')
    with hdf5.archive(path, 'w') as archive:
        archive.create_group('/parameters')
        for name, (mean, error) in VALUES.items():
            alea.ReportedEstimate(count=10, mean=mean, error=error).save(
                archive, RESULTS+'/'+pyalps.hdf5_name_encode(name))
    return path


def test_names_and_components_roundtrip(results_archive):
    with hdf5.archive(results_archive) as archive:
        children = archive.list_children(RESULTS)
        assert sorted(children) == sorted(pyalps.hdf5_name_encode(name) for name in VALUES)
        assert sorted(pyalps.hdf5_name_decode(name) for name in children) == sorted(VALUES)
        encoded_twice = pyalps.hdf5_name_encode('Energy Density &#47; Site')
        assert not archive.is_group(RESULTS+'/'+encoded_twice)
        for stored in children:
            result = alea.read_result(archive, RESULTS+'/'+stored)
            mean, error = VALUES[pyalps.hdf5_name_decode(stored)]
            np.testing.assert_array_equal(result.mean, mean)
            np.testing.assert_array_equal(result.error, error)
            assert result.size == len(mean)
    loaded = {data.props['observable']: data for data in pyalps.loadMeasurements([results_archive])[0]}
    assert loaded.keys() == VALUES.keys()
    for name, (mean, error) in VALUES.items():
        np.testing.assert_array_equal([value.mean for value in loaded[name].y], mean)
        np.testing.assert_array_equal([value.error for value in loaded[name].y], error)


@pytest.mark.parametrize('statistic', ['error', 'count', 'timeseries'])
def test_legacy_measurements_require_explicit_conversion(tmp_path, statistic):
    from pyalps.load import Hdf5Loader
    path = str(tmp_path/'legacy.h5')
    with hdf5.archive(path, 'w') as archive:
        archive.create_group('/parameters')
        archive['/simulation/results/E/mean/value'] = 2.
        field = {'error': 'mean/error', 'count': 'count', 'timeseries': 'timeseries/data'}[statistic]
        archive['/simulation/results/E/'+field] = 1.
    loader = Hdf5Loader()
    with pytest.raises(ValueError, match='alps-hdf5-convert'):
        loader.ReadMeasurementFromFile([path])
    assert loader.h5f.closed


def test_histograms_and_deterministic_values_still_load(tmp_path):
    path = str(tmp_path/'deterministic.h5')
    with hdf5.archive(path, 'w') as archive:
        archive.create_group('/parameters')
        archive['/simulation/results/E/mean/value'] = -2.
        archive['/simulation/results/C/mean/value'] = np.array([1., .5])
        archive['/simulation/results/H/histogram'] = np.array([2, 7, 1])
        archive['/simulation/results/H/@min'] = -.5
        archive['/simulation/results/H/@stepsize'] = .25
    data = {d.props['observable']: d for d in pyalps.loadMeasurements([path])[0]}
    np.testing.assert_array_equal(data['E'].y, [-2.])
    np.testing.assert_array_equal(data['C'].y, [1., .5])
    np.testing.assert_array_equal(data['H'].y, [2, 7, 1])
    np.testing.assert_array_equal(data['H'].x, [-.5, -.25, 0.])
