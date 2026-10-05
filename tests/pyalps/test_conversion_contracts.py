"""Regression checks found by comparing the Boost.Python conversion contracts."""

import os
import subprocess
import sys
import textwrap

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


def test_mcvector_constructor_sizes_errors_before_indexing():
    # Printing/indexing a mean-only vector previously read an empty error
    # vector and crashed the interpreter. Keep native crash checks isolated.
    code = """
        import numpy as np
        from pyalps.alea import MCVectorData
        data = MCVectorData([1.0, 2.0, 3.0])
        np.testing.assert_array_equal(data.error, [0.0, 0.0, 0.0])
        assert data[1].mean == 2.0
        assert data[1].error == 0.0
        assert '2' in repr(data)
        assert '2.00' in format(data, '.2f')
        try:
            MCVectorData([1.0, 2.0], [0.1])
        except ValueError:
            pass
        else:
            raise AssertionError('mismatched mean/error lengths must be rejected')
    """
    result = subprocess.run(
        [sys.executable, "-X", "faulthandler", "-c", textwrap.dedent(code)],
        capture_output=True, text=True, timeout=30,
        env={**os.environ, "MallocScribble": "1"},
    )
    assert result.returncode == 0, result.stdout + result.stderr


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
