"""Measurement workflows must retain native evidence rather than plot summaries."""
import copy

import numpy as np
import pytest

import pyalps
from pyalps import alea, hdf5
from pyalps.dataset import DataSet

FAMILIES = [prefix + family for prefix in ('', 'Complex')
            for family in ('Mean', 'Variance', 'Covariance', 'Autocorrelation', 'Batch')]
FAMILIES += ['EllipticVariance', 'EllipticCovariance']


def measurement(family, start, stop, name='Energy / site'):
    accumulator = getattr(alea, family + 'Accumulator')(2)
    for i in range(start, stop):
        sample = np.array([i % 7, 2*(i % 7) + i % 3], dtype=float)
        if family.startswith(('Complex', 'Elliptic')):
            sample = sample + 1j*sample[::-1]
        accumulator << sample
    return DataSet.from_result(accumulator.result(), x=np.array([.25, .75]),
                               props={'observable': name, 'L': 2, 'seed': start})


@pytest.mark.parametrize('family', FAMILIES)
def test_merge_and_save_preserve_native_evidence(tmp_path, family):
    first, second = measurement(family, 0, 17), measurement(family, 17, 61)
    inputs = [first, second]
    original = copy.deepcopy(inputs)
    merged = pyalps.mergeDataSets(inputs)
    assert len(inputs) == 2 and inputs[0] is first and inputs[1] is second
    assert merged.props == {'observable': 'Energy / site', 'L': 2}
    assert merged.native_result.count == 61
    np.testing.assert_allclose(merged.native_result.mean,
                              (17*first.native_result.mean + 44*second.native_result.mean)/61)
    for actual, expected in zip(inputs, original):
        assert actual.native_result.count == expected.native_result.count
        np.testing.assert_array_equal(actual.native_result.mean, expected.native_result.mean)
    if family.endswith('Batch'):
        np.testing.assert_array_equal(merged.native_result.batch_counts,
                                      np.concatenate([first.native_result.batch_counts, second.native_result.batch_counts]))
        np.testing.assert_array_equal(merged.native_result.batch_sums,
                                      np.concatenate([first.native_result.batch_sums, second.native_result.batch_sums]))
    if family.endswith('Covariance'):
        sample = np.array([[i % 7, 2*(i % 7) + i % 3] for i in range(61)], dtype=float)
        if family == 'Covariance':
            np.testing.assert_allclose(merged.native_result.covariance, np.cov(sample, rowvar=False))
        elif family == 'ComplexCovariance':
            sample = sample + 1j*sample[:, ::-1]
            # ALEA uses E[(x-mean)(x-mean)^H].
            np.testing.assert_allclose(merged.native_result.covariance, np.cov(sample, rowvar=False))
    path = tmp_path/'measurements.h5'
    with hdf5.archive(path, 'w') as ar:
        ar['/unrelated'] = 7
    pyalps.saveMeasurements([merged], path)
    with hdf5.archive(path) as ar:
        assert ar['/unrelated'] == 7
    loaded = pyalps.loadMeasurements([str(path)])[0][0]
    assert loaded.props['observable'] == 'Energy / site'
    np.testing.assert_array_equal(loaded.x, merged.x)
    actual, expected = loaded.native_result, merged.native_result
    assert type(actual) is type(expected) and actual.count == expected.count
    for field in ('mean', 'error', 'variance', 'covariance', 'batch_sums', 'batch_counts', 'count2'):
        if hasattr(expected, field):
            np.testing.assert_array_equal(getattr(actual, field), getattr(expected, field))
    if family.endswith('Autocorrelation'):
        assert actual.levels == expected.levels
        for i in range(actual.levels):
            np.testing.assert_array_equal(actual.level(i).variance, expected.level(i).variance)


def test_group_by_name_and_coordinate_checks():
    left = [measurement('Batch', 0, 20, name) for name in ('E', 'M')]
    right = [measurement('Batch', 20, 50, name) for name in ('M', 'E')]
    merged = pyalps.mergeMeasurements([left, right])
    assert {d.props['observable'] for d in merged} == {'E', 'M'}
    assert all(d.native_result.count == 50 for d in merged)
    with pytest.raises(ValueError, match='empty'):
        pyalps.mergeDataSets([])
    right[1].x = np.array([1., 2.])
    with pytest.raises(ValueError, match='x values'):
        pyalps.mergeDataSets([left[0], right[1]])
    with pytest.raises(ValueError, match='families'):
        pyalps.mergeDataSets([left[0], measurement('Covariance', 0, 20)])


def test_invalid_save_is_preflighted_and_plot_edits_are_not_silently_lost(tmp_path):
    valid = measurement('Batch', 0, 30)
    edited = copy.deepcopy(valid)
    edited.y[0] *= 2
    for invalid in (edited, DataSet([0], [1.], {'observable': 'summary'})):
        target = tmp_path/'must-not-be-created.h5'
        with pytest.raises(ValueError):
            pyalps.saveMeasurements([valid, invalid], target)
        assert not target.exists()
        with pytest.raises(ValueError):
            pyalps.mergeDataSets([valid, invalid])
    with pytest.raises(ValueError, match='duplicate'):
        pyalps.saveMeasurements([valid, valid], tmp_path/'duplicate.h5')
    assert not (tmp_path/'duplicate.h5').exists()


def test_merge_from_files_uses_native_result_path(tmp_path):
    files = []
    for start, stop in ((0, 19), (19, 55)):
        path = tmp_path / f"run-{start}.h5"
        pyalps.saveMeasurements([measurement('Batch', start, stop)], path)
        files.append(str(path))
    merged = pyalps.mergeMeasurementsFromFiles(files)
    assert len(merged) == 1 and merged[0].native_result.count == 55
    np.testing.assert_array_equal(merged[0].x, [.25, .75])


@pytest.mark.parametrize('coordinates', [np.array(['--', 'spin']), np.array([[0, 1], [1, 0]])])
def test_component_labels_roundtrip_without_evaluation(tmp_path, coordinates):
    data = measurement('Batch', 0, 30)
    data.x = coordinates
    target = tmp_path/'labels.h5'
    pyalps.saveMeasurements([data], target)
    loaded = pyalps.loadMeasurements([str(target)])[0][0]
    np.testing.assert_array_equal(loaded.x, coordinates)


def test_labels_are_literals_and_legacy_site_pairs_still_work(tmp_path):
    from pyalps.load import parse_labels
    np.testing.assert_array_equal(parse_labels(['0--1', '0--2']), [1, 2])
    np.testing.assert_array_equal(parse_labels(['0--1', '1--2']), [[0, 1], [1, 2]])
    marker = tmp_path/'must-not-exist'
    with pytest.raises((ValueError, SyntaxError)):
        parse_labels([f"__import__('pathlib').Path({str(marker)!r}).touch()"])
    assert not marker.exists()
