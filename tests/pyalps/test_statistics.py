"""Modern ALEA uses the native scientific core and its canonical checkpoint codec."""
import copy

import h5py
import numpy as np
import pytest

from pyalps import alea, hdf5, ngs
import pyalps


@pytest.mark.parametrize("complex_values", [False, True])
@pytest.mark.parametrize("components", [1, 3])
@pytest.mark.parametrize("cut", [11, 35, 128])
def test_batch_statistics_and_mid_batch_continuation(tmp_path, complex_values, components, cut):
    accumulator_type = alea.ComplexBatchAccumulator if complex_values else alea.BatchAccumulator
    result_type = alea.ComplexBatchResult if complex_values else alea.BatchResult
    x = np.arange(257) % 17 - 8.0
    samples = np.column_stack((x, 2*x, np.sin(x)))[:, :components]
    if complex_values:
        samples = samples + 1j * np.column_stack((np.cos(x), -x, x % 3))[:, :components]
    full = accumulator_type(components, num_batches=8, base_size=3)
    stopped = accumulator_type(components, num_batches=8, base_size=3)
    for i, sample in enumerate(samples):
        full << (sample[0] if components == 1 else sample)
        if i < cut:
            stopped << sample
    filename = tmp_path / "checkpoint.h5"
    hdf5.save_checkpoint(filename, lambda archive: stopped.save(archive, "/measurements"))
    with hdf5.archive(filename) as archive:
        resumed = accumulator_type.read(archive, "/measurements")
    for sample in samples[cut:]:
        resumed << sample
    expected, actual = full.result(), resumed.result()
    assert actual.count == len(samples)
    for field in ("mean", "error", "covariance", "batch_sums", "batch_counts"):
        np.testing.assert_array_equal(getattr(actual, field), getattr(expected, field))
    np.testing.assert_allclose(actual.mean, samples.mean(axis=0), atol=1e-14)
    # Check the weighted estimator independently, including complex covariance.
    counts = actual.batch_counts.astype(float)
    residual = actual.batch_sums[counts > 0] / counts[counts > 0, None] - actual.mean
    counts = counts[counts > 0]
    count2 = counts @ counts
    covariance = residual.T @ (counts[:, None] * residual.conj()) / (counts.sum()-count2/counts.sum())
    np.testing.assert_allclose(actual.covariance, covariance, rtol=1e-12, atol=1e-14)
    np.testing.assert_allclose(actual.error, np.sqrt(covariance.diagonal().real*count2/counts.sum()**2),
                               rtol=1e-12, atol=1e-14)
    assert actual.count2 == count2
    assert actual.observations == pytest.approx(counts.sum()**2/count2)
    with hdf5.archive(filename, "a") as archive:
        archive["/result"] = actual
        loaded = result_type.read(archive, "/result")
        np.testing.assert_array_equal(loaded.batch_sums, actual.batch_sums)
        np.testing.assert_array_equal(loaded.covariance, actual.covariance)
    with h5py.File(filename) as archive:
        assert archive["result"].attrs["kind"] == 5
        assert archive["measurements"].attrs["kind"] == 6
        assert archive["result/batch/sum"].shape == actual.batch_sums.shape
        assert archive["result/batch/sum"].dtype == np.dtype("complex128" if complex_values else "float64")
        assert archive["result/batch/count"].dtype == np.dtype("uint64")


@pytest.mark.parametrize("complex_values", [False, True])
@pytest.mark.parametrize("components", [1, 3])
def test_analysis_loader_preserves_modern_means_and_errors(tmp_path, complex_values, components):
    accumulator_type = alea.ComplexBatchAccumulator if complex_values else alea.BatchAccumulator
    accumulator = accumulator_type(components, num_batches=8, base_size=3)
    for i in range(37):
        sample = np.arange(components, dtype=float) + i % 11
        if complex_values:
            sample = sample + 1j * (i % 7 - np.arange(components))
        accumulator << sample
    result = accumulator.result()
    filename = str(tmp_path / "modern-results.h5")
    name = "Energy / per site"
    path = "/simulation/results/" + pyalps.hdf5_name_encode(name)
    with hdf5.archive(filename, "w") as archive:
        archive["/parameters"] = ngs.params({"L": components})
        result.save(archive, path)
        accumulator_type(1).result().save(archive, "/simulation/results/empty")
    loaded = pyalps.loadMeasurements([filename])
    assert len(loaded) == 1 and len(loaded[0]) == 1
    measured = loaded[0][0]
    assert measured.props["observable"] == name
    np.testing.assert_array_equal([value.mean for value in measured.y], result.mean)
    np.testing.assert_array_equal([value.error for value in measured.y], result.error)
    np.testing.assert_array_equal(measured.x, np.arange(components))


@pytest.mark.parametrize("histogram", [False, True])
def test_ctint_publishes_native_statistics_for_python_analysis(tmp_path, monkeypatch, histogram):
    from pyalps import ctint

    monkeypatch.chdir(tmp_path)
    filename = str(tmp_path / "ctint-results.h5")
    run = ctint.prepare(
        {"BETA": 2., "U": 0., "MU": 0., "ALPHA": -.01, "N": 8,
         "NMATSUBARA": 4, "SWEEPS": 37, "THERMALIZATION": 2,
         "MEASUREMENT_PERIOD": 1, "HISTOGRAM_MEASUREMENT": histogram},
        input={"atomic": True}, output={"results": filename}, execution={"bins": 8})
    ctint.solve(run)
    measured = {entry.props["observable"]: entry for entry in pyalps.loadMeasurements([filename])[0]}
    with hdf5.archive(filename) as archive:
        assert archive["/run_config/application"] == "ctint"
        for name, entry in measured.items():
            path = "/simulation/results/" + pyalps.hdf5_name_encode(name)
            assert archive[path + "/@version"] == 1 and archive[path + "/@kind"] == 5
            result = alea.BatchResult.read(archive, path)
            assert result.count > 1
            np.testing.assert_array_equal([value.mean for value in entry.y], result.mean)
            np.testing.assert_array_equal([value.error for value in entry.y], result.error)
        np.testing.assert_allclose(archive["/G_omega/0/mean/value"],
                                   -2j / ((2 * np.arange(4) + 1) * np.pi), atol=1e-12)
    assert measured["Sign"].y[0].mean == 1.
    assert measured["Sign"].y[0].error == 0.
    # Fermion occupation is idempotent: n_sigma**2 = n_sigma. With the
    # solver's Sz = n_up - n_down convention, the atomic spin square is 1/2.
    if histogram:
        assert measured["Sz_0"].y[0].mean == 0.
        assert measured["Sz2_0"].y[0].mean == pytest.approx(.5)
        assert measured["Sz0_Sz0"].y[0].mean == pytest.approx(.5)
    else:
        np.testing.assert_allclose([value.mean for value in measured["n_i n_j"].y], [.5, .25, .25, .5])
    assert {path.name for path in tmp_path.iterdir()} == {"ctint-results.h5"}


def test_cthyb_publishes_batched_and_component_results_for_python_analysis(tmp_path, monkeypatch):
    from pyalps import cthyb

    monkeypatch.chdir(tmp_path)
    delta = tmp_path / "delta.dat"
    delta.write_text("".join(f"{i} -0.5 -0.5\n" for i in range(9)))
    filename = str(tmp_path / "cthyb-results.h5")
    run = cthyb.prepare(
        {"BETA": 2., "U": 0., "MU": 0., "N_ORBITALS": 2, "N_TAU": 8,
         "N_MEAS": 3, "THERMALIZATION": 0, "SWEEPS": 37,
         "N_MATSUBARA": 4, "N_LEGENDRE": 4, "N_nn": 4, "N_w2": 2, "N_W": 1,
         "MEASURE_freq": True, "MEASURE_legendre": True, "MEASURE_nn": True,
         "MEASURE_nnt": True, "MEASURE_nnw": True, "MEASURE_g2w": True,
         "MEASURE_h2w": True, "MEASURE_sector_statistics": True},
        input={"delta": str(delta)}, output={"results": filename, "base_path": "/pilot"},
        execution={"bins": 8})
    cthyb.solve(run)
    measured = {entry.props["observable"]: entry for entry in
                pyalps.loadMeasurements([filename], respath="/pilot/simulation/results")[0]}
    kinds = set()
    with hdf5.archive(filename) as archive:
        assert archive["/run_config/application"] == "cthyb"
        for name, entry in measured.items():
            path = "/pilot/simulation/results/" + pyalps.hdf5_name_encode(name)
            kind = archive[path + "/@kind"]
            kinds.add(kind)
            assert archive[path + "/@version"] == 1
            result_type = alea.VarianceResult if kind == 2 else alea.BatchResult
            result = result_type.read(archive, path)
            assert result.count == 37
            np.testing.assert_array_equal([value.mean for value in entry.y], result.mean)
            np.testing.assert_array_equal([value.error for value in entry.y], result.error)
            assert np.isfinite(result.mean).all() and np.isfinite(result.error).all()
            if name.startswith(("g2w_", "h2w_")):
                assert kind == 2 and not archive.is_group(path + "/batch")
                assert result.variance.shape == result.mean.shape
                assert result.count2 == result.count
        for orbital in range(2):
            g = archive[f"/pilot/G_tau/{orbital}/mean/value"]
            np.testing.assert_allclose(g[0] + g[-1], -1., atol=1e-14)
            for observable in ("G_tau", "F_tau"):
                path = f"/pilot/{observable}/{orbital}/mean/"
                error = archive[path + "error"]
                covariance = archive[path + "covariance"].reshape(9, 9)
                np.testing.assert_allclose(covariance.diagonal(), error**2, atol=1e-14)
                if observable == "G_tau":
                    np.testing.assert_allclose(covariance[0], -covariance[-1], atol=1e-14)
                else:
                    # U=0 makes F identically zero; copying G's covariance is incorrect.
                    np.testing.assert_array_equal(archive[path + "value"], np.zeros(9))
                    np.testing.assert_array_equal(covariance, np.zeros((9, 9)))
    assert kinds == {2, 5}
    assert measured["Sign"].y[0].mean == 1.
    assert measured["Sign"].y[0].error == 0.
    assert {path.name for path in tmp_path.iterdir()} == {"delta.dat", "cthyb-results.h5"}


def test_batch_rejected_sample_and_failed_load_preserve_state(tmp_path):
    accumulator = alea.BatchAccumulator(2, num_batches=8)
    for i in range(8):
        accumulator << np.array([i, 2*i], dtype=float)
    before = accumulator.result()
    for sample in (np.ones(3), np.ones((1, 2)), np.array([1j, 2j]), ["1", "2"]):
        with pytest.raises((RuntimeError, ValueError, TypeError)):
            accumulator << sample
        np.testing.assert_array_equal(accumulator.result().batch_sums, before.batch_sums)
        np.testing.assert_array_equal(accumulator.result().batch_counts, before.batch_counts)
    filename = tmp_path / "checkpoint.h5"
    hdf5.save_checkpoint(filename, lambda archive: accumulator.save(archive))
    with h5py.File(filename, "a") as archive:
        archive["cursor/level"][()] = np.uint64(64)
    with hdf5.archive(filename) as archive:
        archive.set_context("/caller")
        with pytest.raises(RuntimeError):
            accumulator.load(archive, "/")
        assert archive.context == "/caller"
    np.testing.assert_array_equal(accumulator.result().batch_sums, before.batch_sums)
    np.testing.assert_array_equal(accumulator.result().batch_counts, before.batch_counts)
    # Returned arrays belong to the result, not to a borrowed accumulator store.
    sums = before.batch_sums
    sums.fill(99)
    np.testing.assert_array_equal(accumulator.result().batch_sums, before.batch_sums)


def test_python_checkpoint_publication_closes_callbacks_and_rolls_back(tmp_path):
    filename = tmp_path / "state.h5"
    retained = []

    def save(archive):
        retained.append(archive)
        archive["/state"] = 1

    hdf5.save_checkpoint(filename, save)
    assert not retained[0].is_open
    with pytest.raises(hdf5.ArchiveClosed):
        retained[0]["/state"]
    before = filename.read_bytes()

    def fail(archive):
        archive["/state"] = 2
        raise ValueError("failed scientific save")

    with pytest.raises(ValueError, match="failed scientific save"):
        hdf5.save_checkpoint(filename, fail)
    assert filename.read_bytes() == before
    with pytest.raises(ValueError):
        hdf5.save_checkpoint(tmp_path / "missing.h5", fail)
    assert set(tmp_path.iterdir()) == {filename}


@pytest.mark.parametrize("name", ["mt19937", "lagged_fibonacci607"])
@pytest.mark.parametrize("state", ["invalid engine", "trailing"])
def test_random_checkpoint_failure_preserves_the_stream(tmp_path, state, name):
    random = ngs.random01(17, name)
    for _ in range(31):
        random()
    expected = copy.deepcopy(random)
    filename = tmp_path / "random.h5"
    with hdf5.archive(filename, "w") as archive:
        archive["/"] = random
        archive["/engine"] = archive["/engine"] + " junk" if state == "trailing" else state
        with pytest.raises(RuntimeError, match="invalid random01 checkpoint"):
            random.load(archive)
    assert [random() for _ in range(64)] == [expected() for _ in range(64)]


@pytest.mark.parametrize("name", ["mt19937", "lagged_fibonacci607"])
def test_random_engines_checkpoint_and_copy(tmp_path, name):
    random = ngs.random01(37, name)
    values = [random() for _ in range(1001)]
    assert all(0 <= value < 1 for value in values)
    expected = copy.deepcopy(random)
    with hdf5.archive(tmp_path / "random.h5", "w") as archive:
        random.save(archive)
        restored = ngs.random01(0)
        restored.load(archive)
    assert restored.name == name
    assert [restored() for _ in range(1001)] == [expected() for _ in range(1001)]
