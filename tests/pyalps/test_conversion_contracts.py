"""Regression checks found by comparing the Boost.Python conversion contracts."""


import numpy as np
import pytest

from pyalps import alea




@pytest.mark.parametrize("layout", ["readonly", "strided", "fortran"])
def test_array_consumers_do_not_require_writable_samples(layout):
    values = np.arange(24.0).reshape(4, 6)
    if layout == "readonly":
        values.flags.writeable = False
    elif layout == "strided":
        values = values[:, ::2]
    else:
        values = np.asfortranarray(values)
    np.testing.assert_allclose(alea.mean(values[0]), values[0].mean())
    np.testing.assert_allclose(alea.mean(values), values.mean(axis=0))
    native = alea.BatchAccumulator(values.shape[1])
    for row in values:
        native << row
    np.testing.assert_allclose(native.result().mean, values.mean(axis=0))


@pytest.mark.parametrize('vector', [False, True])
def test_full_timeseries_uses_numpy_storage_and_native_statistics(tmp_path, vector):
    from pyalps import hdf5
    samples = np.arange(257, dtype=float) % 19 - 9
    if vector:
        samples = np.column_stack((samples, samples*2+1))
    statistics = alea.BatchAccumulator(2 if vector else 1, num_batches=8)
    diagnostics = alea.AutocorrelationAccumulator(2 if vector else 1)
    for sample in samples:
        statistics << sample
        diagnostics << sample
    filename = tmp_path/'samples.h5'
    def save(ar):
        ar['/samples'] = samples
        statistics.result().save(ar, '/result')
        diagnostics.result().save(ar, '/diagnostics')
    hdf5.save_checkpoint(filename, save)
    with hdf5.archive(filename) as ar:
        np.testing.assert_array_equal(ar['/samples'], samples)
        result = alea.read_result(ar, '/result')
        analysis = alea.read_result(ar, '/diagnostics')
    np.testing.assert_allclose(result.mean, np.atleast_1d(samples.mean(axis=0)))
    assert result.count == analysis.count == len(samples)
