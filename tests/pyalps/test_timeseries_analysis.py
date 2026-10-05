"""Numerical contracts for chronological NumPy sample analysis."""
import numpy as np
import pytest
from pyalps import alea


@pytest.mark.parametrize('vector', [False, True])
def test_correlation_matches_direct_lag_products(vector):
    x = np.random.default_rng(19).normal(size=(127, 3) if vector else 127)
    centered = x - x.mean(axis=0)
    expected = np.asarray([(centered[:-k]*centered[k:]).sum(axis=0) /
                           ((len(x)-k)*x.var(axis=0, ddof=1)) for k in range(1, len(x))])
    np.testing.assert_allclose(alea.autocorrelation(x, _distance=-1), expected, atol=1e-14)
    np.testing.assert_allclose(alea.error(x), np.sqrt(x.var(axis=0, ddof=1)/len(x)))
    native = alea.AutocorrelationAccumulator(3 if vector else 1)
    for sample in x:
        native << sample
    np.testing.assert_allclose(alea.error(x, 'binning'), native.result().error if vector else native.result().error[0])
    np.testing.assert_allclose(alea.running_mean(x), [x[:i].mean(axis=0) for i in range(1, len(x)+1)])
    np.testing.assert_allclose(alea.reverse_running_mean(x), [x[i:].mean(axis=0) for i in range(len(x))])


def test_exponential_fit_uses_absolute_lags_and_integrated_tail():
    lags = np.arange(1, 101)
    corr = 1.4*np.exp(-0.12*lags)
    for bounds in ({'_from': 7, '_to': 25}, {'_max': .8, '_min': .2}, {'_from': -94, '_to': -75}):
        fit = alea.exponential_autocorrelation_time(corr, **bounds)
        np.testing.assert_allclose(fit, [1.4, -.12], rtol=1e-13)
        expected = corr.sum() + 1.4/.12*np.exp(-.12*100.5)
        assert alea.integrated_autocorrelation_time(corr, fit) == pytest.approx(expected)
    assert alea.integrated_autocorrelation_time(corr) == pytest.approx(corr.sum())


def test_cut_views_thresholds_and_dataset():
    x = np.array([10., 8., 5., 3., 1.])
    np.testing.assert_array_equal(alea.cut_head(x, _limit=.5), [3., 1.])
    np.testing.assert_array_equal(alea.cut_tail(x, _limit=.5), [10., 8., 5.])
    for result in (alea.cut_head(x, _distance=2), alea.cut_tail(x, _distance=2)):
        assert np.shares_memory(result, x)
    assert len(alea.cut_head(x, _distance=len(x))) == 0
    assert len(alea.cut_tail(x, _distance=len(x))) == 0
    np.testing.assert_array_equal(alea.cut_tail(x, _distance=0), x)
    dataset = alea.make_dataset(x)
    np.testing.assert_array_equal(dataset.x, [1, 2, 3, 4, 5])
    np.testing.assert_array_equal(dataset.y, x)


def test_limit_includes_first_crossing_in_any_component():
    x = np.random.default_rng(72).normal(size=(100, 2)).cumsum(axis=0)
    full = alea.autocorrelation(x, _distance=99)
    first = np.flatnonzero((full < .5*full[0]).any(axis=1))[0]
    np.testing.assert_allclose(alea.autocorrelation(x, _limit=.5), full[:first+1])


@pytest.mark.parametrize('x', [[], [float('nan')], [float('inf')], [1j], [['x']], np.zeros((2,0)), np.zeros((2,2,2))])
def test_invalid_mean_samples(x):
    with pytest.raises(ValueError):
        alea.mean(x)


@pytest.mark.parametrize('call', [
    lambda: alea.autocorrelation([1., 1.], _distance=1),
    lambda: alea.autocorrelation([1., 2.]),
    lambda: alea.autocorrelation([1., 2.], _distance=1, _limit=.5),
    lambda: alea.autocorrelation([1., 2.], _distance=2),
    lambda: alea.variance([1.]),
    lambda: alea.error([1., 2.], 'unknown'),
    lambda: alea.cut_head([1., 2.], _distance=3),
    lambda: alea.cut_tail([[1., 2.]], _limit=.5),
    lambda: alea.exponential_autocorrelation_time([1., 0.], _from=1, _to=2),
    lambda: alea.exponential_autocorrelation_time([1., .5], _from=1, _to=1),
    lambda: alea.integrated_autocorrelation_time([1., .5], (1., 0.)),
])
def test_invalid_analysis_requests(call):
    with pytest.raises(ValueError):
        call()
