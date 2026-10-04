"""Modern ALEA uses the native scientific core and its canonical checkpoint codec."""
import copy

import h5py
import numpy as np
import pytest

from pyalps import alea, hdf5, ngs


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


@pytest.mark.parametrize("state", ["invalid engine", "trailing"])
def test_random_checkpoint_failure_preserves_the_stream(tmp_path, state):
    random = ngs.random01(17)
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
