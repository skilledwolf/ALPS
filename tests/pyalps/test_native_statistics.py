"""Native statistical evidence, continuation and analysis through Python."""
import copy

import h5py
import numpy as np
import pytest

import pyalps
from pyalps import alea, hdf5


FAMILIES = [prefix + family for prefix in ("", "Complex")
            for family in ("Mean", "Variance", "Covariance", "Autocorrelation", "Batch")]
FAMILIES += ["EllipticVariance", "EllipticCovariance"]


@pytest.mark.parametrize("family", FAMILIES)
def test_native_estimator_checkpoint_and_loader(tmp_path, family):
    cls = getattr(alea, family + "Accumulator")
    kwargs = (dict(num_batches=8, base_size=3) if family.endswith("Batch") else
              {} if family.endswith("Mean") else dict(batch_size=3))
    acc = cls(2, **kwargs)
    values = np.column_stack((np.arange(281) % 17 - 8., np.arange(281) % 11 - 5.))
    if family.startswith(("Complex", "Elliptic")):
        values = values + 1j * (values[:, ::-1] + 2.)
    for x in values[:37]:
        acc << x
    assert acc.count == 37
    filename = str(tmp_path / "result.h5")
    with hdf5.archive(filename, "w") as ar:
        acc.save(ar, "/checkpoint")
    with hdf5.archive(filename, "r") as ar:
        restored = cls.read(ar, "/checkpoint")
    for x in values[37:]:
        acc << x
        restored << x
    expected, actual = acc.result(), restored.result()
    assert actual.count == len(values)
    for field in ("mean", "error", "variance", "covariance", "count2", "batch_sums", "batch_counts"):
        if hasattr(expected, field):
            np.testing.assert_array_equal(getattr(actual, field), getattr(expected, field))
    np.testing.assert_allclose(actual.mean, values.mean(axis=0), atol=1e-14)
    with hdf5.archive(filename, "a") as ar:
        actual.save(ar, "/simulation/results/value")
    with h5py.File(filename, "a") as ar:
        ar.require_group("parameters")
    data = pyalps.loadMeasurements([filename])[0][0]
    assert type(data.native_result) is type(actual)
    np.testing.assert_array_equal(data.native_result.mean, actual.mean)
    if hasattr(actual, "covariance"):
        np.testing.assert_array_equal(data.native_result.covariance, actual.covariance)
    for d in (data, copy.deepcopy(data)):
        assert d.native_result.count == actual.count


@pytest.mark.parametrize("family", FAMILIES)
def test_independent_merge_and_join_keep_weights(family):
    cls = getattr(alea, family + "Accumulator")
    left, right = cls(1), cls(1)
    for i in range(10):
        left << float(i)
    for i in range(10, 30):
        right << float(i)
    a, b = left.result(), right.result()
    pooled = alea.merge([a, b])
    assert pooled.count == 30
    np.testing.assert_allclose(pooled.mean, [14.5])
    joined = a.join(a)
    np.testing.assert_allclose(joined.mean, [4.5, 4.5])
    if hasattr(joined, "count2"):
        assert joined.count2 == a.count2
    if hasattr(pooled, "variance"):
        expected = np.var(np.arange(30.), ddof=1)
        if family.startswith("Elliptic"):
            np.testing.assert_allclose(pooled.variance[:, 0, 0], [expected])
            np.testing.assert_array_equal(pooled.variance[:, 1, 1], [0.])
        else:
            np.testing.assert_allclose(pooled.variance, [expected])
    assert a.count == 10 and b.count == 20
    with pytest.raises(ValueError):
        alea.merge([])


def test_native_covariance_propagation_and_mean_test():
    acc = alea.CovarianceAccumulator(2)
    batch = alea.BatchAccumulator(2, num_batches=16)
    for x in np.arange(1., 257.):
        acc << [x, 2 * x]
        batch << [x, 2 * x]
    result = acc.result()
    difference = result.transform(lambda x: np.array([x[1] - 2*x[0]]))
    np.testing.assert_allclose(difference.mean, [0.], atol=1e-12)
    np.testing.assert_allclose(difference.error, [0.], atol=1e-6)
    ratio = batch.result().transform(lambda x: np.array([x[1]/x[0]]))
    np.testing.assert_allclose(ratio.mean, [2.])
    np.testing.assert_allclose(ratio.error, [0.], atol=1e-12)
    assert result.test_mean(result.mean).pvalue == 1.
    assert result.test_mean(result).pvalue == 1.
    with pytest.raises(ValueError):
        result.transform(lambda x: x, output_size=1)
    with pytest.raises(ValueError):
        result.transform(lambda x: np.array([x[0]]), method="jackknife")


def test_autocorrelation_hierarchy_and_elliptic_uncertainty():
    acc = alea.AutocorrelationAccumulator()
    for x in np.repeat(np.arange(64) % 7 - 3., 128):
        acc << x
    result = acc.result()
    assert result.tau_available
    assert result.levels > 1
    assert result.tau[0] > 0
    assert result.level(0).count == result.count
    assert result.error[0] > result.level(0).error[0]
    with pytest.raises(IndexError):
        result.level(result.levels)
    joint = alea.EllipticVarianceAccumulator()
    for _ in range(10000):
        joint << .64 + 1j
    assert joint.result().variance.shape == (1, 2, 2)
    np.testing.assert_array_equal(alea.ratio_real_imag(joint.result()).error, [0.])


def test_complex_real_components_retain_nonlinear_and_elliptic_evidence():
    batch = alea.ComplexBatchAccumulator(1, num_batches=8)
    cov = alea.EllipticCovarianceAccumulator(1)
    for x in np.arange(1., 65.):
        batch << x + 2j*x
        cov << x + 2j*x
    for result in (batch.result().real_components(), cov.result().real_components()):
        np.testing.assert_allclose(result.mean, [32.5, 65.])
        transformed = result.transform(lambda x: np.array([x[1]/x[0]]))
        np.testing.assert_allclose(transformed.mean, [2.])
        np.testing.assert_allclose(transformed.error, [0.], atol=1e-7)
    circular = alea.ComplexCovarianceAccumulator()
    circular << 1j
    with pytest.raises(ValueError, match="elliptic"):
        circular.result().real_components()


def test_mcbase_preserves_heterogeneous_native_estimators(tmp_path):
    from pyalps import ngs

    class Simulation(ngs.mcbase):
        def __init__(self):
            super().__init__({"SEED": 123}, 0)

        def update(self):
            pass

        def measure(self):
            pass

        def fraction_completed(self):
            return 0.

    uninterrupted = Simulation()
    for family in FAMILIES:
        uninterrupted.measurements[family] = getattr(alea, family + "Accumulator")(2)
    samples = np.column_stack((np.arange(47) % 11 - 5., np.arange(47) % 7 - 3.))

    def sample(sim, values):
        for x in values:
            for name in sim.measurements:
                sim.measurements[name] << (x + 1j * x[::-1] if name.startswith(("Complex", "Elliptic")) else x)

    sample(uninterrupted, samples[:19])
    filename = str(tmp_path / "mixed.h5")
    with hdf5.archive(filename, "w") as ar:
        uninterrupted.save(ar)
    resumed = Simulation()
    with hdf5.archive(filename, "r") as ar:
        resumed.load(ar)
    sample(uninterrupted, samples[19:])
    sample(resumed, samples[19:])
    expected, actual = uninterrupted.collectResults(), resumed.collectResults()
    assert expected.keys() == actual.keys()
    for name in FAMILIES:
        assert type(actual[name]) is getattr(alea, name + "Result")
        for field in ("mean", "error", "covariance", "variance", "batch_sums", "batch_counts"):
            if hasattr(expected[name], field):
                np.testing.assert_array_equal(getattr(expected[name], field), getattr(actual[name], field))
    # Type mismatches reject before replacing any registered state.
    resumed.measurements["Mean"] = alea.CovarianceAccumulator(2)
    retained = resumed.measurements["Mean"]
    with hdf5.archive(filename, "r") as ar:
        with pytest.raises(Exception):
            resumed.load(ar)
    assert resumed.measurements["Mean"] is retained
    assert resumed.collectResults()["Batch"].count == 47
    with hdf5.archive(str(tmp_path / "results.h5"), "w") as ar:
        ngs.saveResults(expected, uninterrupted.parameters, ar, "/simulation/results")


@pytest.mark.parametrize('complex_values', [False, True])
def test_error_convergence_preserves_rising_error_evidence(tmp_path, complex_values):
    cls = alea.ComplexAutocorrelationAccumulator if complex_values else alea.AutocorrelationAccumulator
    acc = cls(3)
    for i in range(8192):
        acc << np.full(3, (-1.)**i, dtype=complex if complex_values else float)
    filename = tmp_path/'plateau.h5'
    with hdf5.archive(filename, 'w') as ar:
        acc.result().save(ar, '/result')
    # Construct valid moment summaries with prescribed errors. All four
    # compared levels have >=1024 batches. A later plateau must not erase the
    # earlier uncertain/not-converged evidence in components 1 and 2.
    errors = [np.array([1., .85, .7]), np.ones(3), np.ones(3), np.ones(3)]
    with h5py.File(filename, 'a') as ar:
        for level, error in enumerate(errors):
            ar[f'result/level/{level}/var'][...] = error**2 * (8192 // 2**level)
    with hdf5.archive(filename) as ar:
        result = alea.read_result(ar, '/result')
    assert result.converged_errors.dtype.kind == 'i'
    np.testing.assert_array_equal(result.converged_errors, [0, 1, 2])


def test_error_convergence_requires_four_populated_levels():
    acc = alea.AutocorrelationAccumulator(2)
    np.testing.assert_array_equal(acc.result().converged_errors, [1, 1])
    for i in range(8191):
        acc << [1., -1. if i < 4096 else 1.]
    np.testing.assert_array_equal(acc.result().converged_errors, [1, 1])
    acc << [1., 1.]
    np.testing.assert_array_equal(acc.result().converged_errors, [0, 2])
    acc << [np.nan, 1.]
    assert acc.result().converged_errors[0] == 1
