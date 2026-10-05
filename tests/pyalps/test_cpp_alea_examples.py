"""The C++ and Python analysis lessons consume the same raw chronological data."""
import os
from pathlib import Path
import re
import subprocess
import sys

import numpy as np
import pytest
from pyalps import alea, hdf5


def test_cpp_analysis_examples(tmp_path):
    directory = os.environ.get('ALPS_ALEA_EXAMPLES_DIR')
    if not directory:
        pytest.skip('Set ALPS_ALEA_EXAMPLES_DIR to the built C++ ALEA examples')
    directory = Path(directory).resolve()
    generator = Path(__file__).resolve().parents[2]/'tutorials/00-examples/alea/generate_samples.py'
    subprocess.run([sys.executable, str(generator)], cwd=tmp_path, check=True)
    with hdf5.archive(str(tmp_path/'timeseries.h5')) as archive:
        energy, magnetization = archive['/samples/E'], archive['/samples/m']
    output = {}
    for example in ('mean', 'variance', 'autocorrelation', 'error', 'running_mean'):
        result = subprocess.run([str(directory/f'example_{example}')], cwd=tmp_path,
                                check=True, capture_output=True, text=True)
        output[example] = result.stdout
    correlations = alea.autocorrelation(magnetization, _distance=len(magnetization)-1)
    native_correlations = np.fromstring(output['autocorrelation'].splitlines()[1], sep=' ')
    np.testing.assert_allclose(native_correlations, correlations, atol=1e-11, rtol=1e-10)
    # Independently check normalization and lag indexing with direct dot products.
    centered = magnetization-magnetization.mean()
    for lag in (1, 7, 105, len(magnetization)-1):
        expected = np.dot(centered[:-lag], centered[lag:])/((len(centered)-lag)*np.var(centered, ddof=1))
        np.testing.assert_allclose(native_correlations[lag-1], expected, atol=1e-11)
    fit = alea.exponential_autocorrelation_time(correlations, _max=.8, _min=.2)
    actual_fit = re.search(r'fit is: (\S+) \* exp\((\S+) \* t\)', output['autocorrelation'])
    np.testing.assert_allclose([float(x) for x in actual_fit.groups()], fit, rtol=1e-10)
    tau = alea.integrated_autocorrelation_time(alea.cut_tail(correlations, _limit=.2), fit)
    estimates = dict(line.split(': ') for line in output['error'].splitlines())
    independent = np.std(magnetization, ddof=1)/np.sqrt(len(magnetization))
    expected_error = independent*np.sqrt(1+2*tau)
    np.testing.assert_allclose(float(estimates['uncorrelated']), independent, rtol=1e-12)
    np.testing.assert_allclose(float(estimates['with correlation time']), expected_error, rtol=1e-10)
    np.testing.assert_allclose(float(estimates['with binning']), alea.error(magnetization, alea.binning), rtol=1e-12)
    lines = output['running_mean'].splitlines()
    np.testing.assert_allclose(np.fromstring(lines[1], sep=' '), np.cumsum(energy)/np.arange(1, len(energy)+1))
    np.testing.assert_allclose(np.fromstring(lines[3], sep=' '),
                               (np.cumsum(energy[::-1])/np.arange(1, len(energy)+1))[::-1])
    with hdf5.archive(str(tmp_path/'timeseries.h5')) as archive:
        np.testing.assert_array_equal(archive['/samples/E'], energy)
        np.testing.assert_array_equal(archive['/samples/m'], magnetization)
        np.testing.assert_allclose(archive['/analysis/E/mean'], energy.mean(), rtol=1e-12)
        np.testing.assert_allclose(archive['/analysis/E/variance'], energy.var(ddof=1), rtol=1e-12)
        np.testing.assert_allclose(archive['/analysis/m/error'], expected_error, rtol=1e-10)
