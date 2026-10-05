"""Direct ALEA checkpoint replacement and failed-load contracts."""

import h5py
import numpy as np
import pytest

from pyalps import alea, hdf5


def _full_result(tmp_path, vector):
    samples = np.arange(256) % 17 - 8.0
    if vector:
        samples = np.column_stack((samples, 2.0 * samples + 3.0))
    observable = (alea.RealVectorObservable if vector else alea.RealObservable)("energy")
    for sample in samples:
        observable << sample
    filename = str(tmp_path / "alea.h5")
    observable.save(filename)
    value = (alea.MCVectorData if vector else alea.MCScalarData)()
    value.load(filename, "/simulation/results/energy")
    np.testing.assert_allclose(value.mean, samples.mean(axis=0), rtol=1e-12)
    np.testing.assert_allclose(value.variance, samples.var(axis=0, ddof=1), rtol=1e-12)
    # Materialize the optional jackknife cache before the first save.
    assert len(value.jackknife) == len(value.bins) + 1
    value.save(filename, "/result")
    return filename, value


@pytest.mark.parametrize("vector", [False, True])
def test_direct_alea_save_replaces_optional_leaves_only(tmp_path, vector):
    filename, value = _full_result(tmp_path, vector)
    with h5py.File(filename, "a") as archive:
        for parent in ("variance", "tau", "jacknife"):
            archive[f"result/{parent}/application"] = 19
        archive["application/status"] = 23
        archive["result"].attrs["application"] = "retained"

    # Invalidating a previously persisted cache must remove that cache.
    value.discard_bins(0)
    value.save(filename, "/result")
    with h5py.File(filename, "r") as archive:
        assert "result/jacknife/data" not in archive
        assert "result/variance/value" in archive and "result/tau/value" in archive

    replacement = (alea.MCVectorData([11.0, 12.0], [0.25, 0.5]) if vector
                   else alea.MCScalarData(11.0, 0.25))
    replacement.discard_bins(0)
    replacement.save(filename, "/result")
    with h5py.File(filename, "r") as archive:
        for field in ("variance/value", "tau/value", "jacknife/data"):
            assert f"result/{field}" not in archive
        for parent in ("variance", "tau", "jacknife"):
            assert archive[f"result/{parent}/application"][()] == 19
        assert archive["application/status"][()] == 23
        assert archive["result"].attrs["application"] == "retained"
        assert archive["result/count"][()] == 1
        np.testing.assert_array_equal(archive["result/mean/value"][()], replacement.mean)


@pytest.mark.parametrize("vector", [False, True])
def test_direct_alea_load_replaces_absent_statistics_and_history(tmp_path, vector):
    filename, value = _full_result(tmp_path, vector)
    mean = np.array([41.0, 42.0]) if vector else 41.0
    error = np.array([0.5, 0.75]) if vector else 0.5
    with h5py.File(filename, "a") as archive:
        group = archive.create_group("summary")
        group["count"] = np.uint64(7)
        group["mean/value"] = mean
        group["mean/error"] = error
        group.attrs["cannotrebin"] = False
    value.load(filename, "/summary")
    assert value.count == 7
    np.testing.assert_array_equal(value.mean, mean)
    np.testing.assert_array_equal(value.error, error)
    assert value.bins.size == value.jackknife.size == 0
    for field in ("variance", "tau"):
        with pytest.raises(RuntimeError):
            getattr(value, field)
    value.save(filename, "/replaced")
    with h5py.File(filename, "r") as archive:
        assert archive["replaced/timeseries/data"].attrs["binsize"] == 0
        assert archive["replaced/timeseries/data"].attrs["maxbinnum"] == 0
        assert "replaced/variance/value" not in archive
        assert "replaced/tau/value" not in archive


@pytest.mark.parametrize("vector", [False, True])
def test_direct_alea_empty_checkpoint_replaces_populated_state(tmp_path, vector):
    filename, value = _full_result(tmp_path, vector)
    empty = (alea.MCVectorData if vector else alea.MCScalarData)()
    empty.discard_bins(0)
    empty.save(filename, "/result")
    value.load(filename, "/result")
    assert value.count == 0 and value.bins.size == value.jackknife.size == 0
    for field in ("mean", "error", "variance", "tau"):
        with pytest.raises(RuntimeError):
            getattr(value, field)
    with h5py.File(filename, "r") as archive:
        assert archive["result/count"][()] == 0
        for field in ("variance/value", "tau/value", "jacknife/data"):
            assert f"result/{field}" not in archive


@pytest.mark.parametrize("vector", [False, True])
def test_direct_alea_single_raw_observable_has_unavailable_error(tmp_path, vector):
    filename, value = _full_result(tmp_path, vector)
    sample = np.array([0.0, -2.0]) if vector else -2.0
    observable = (alea.RealVectorObservable if vector else alea.RealObservable)("single")
    observable << sample
    observable.save(filename)
    path = "/simulation/results/single"
    with h5py.File(filename, "r") as archive:
        assert archive[f"{path}/count"][()] == 1
        assert f"{path}/mean/error" not in archive
    value.load(filename, path)
    assert value.count == 1 and value.bins.size == 0
    np.testing.assert_array_equal(value.mean, sample)
    np.testing.assert_array_equal(value.error, np.full(np.asarray(sample).shape, np.inf))
    for field in ("variance", "tau"):
        with pytest.raises(RuntimeError):
            getattr(value, field)
    # A result checkpoint always writes its error. A missing canonical result
    # field remains malformed even when its count happens to be one.
    value.save(filename, "/one_result")
    with h5py.File(filename, "a") as archive:
        del archive["one_result/mean/error"]
    with pytest.raises(RuntimeError, match="require an error"):
        value.load(filename, "/one_result")
    np.testing.assert_array_equal(value.mean, sample)
    assert np.all(np.isposinf(value.error))


@pytest.mark.parametrize("vector", [False, True])
@pytest.mark.parametrize("corruption", ["late_type", "missing_error", "zero_binsize",
                                         "jack_count", "component_shape", "malformed_group"])
def test_direct_alea_malformed_load_preserves_complete_state(tmp_path, vector, corruption):
    filename, value = _full_result(tmp_path, vector)
    fields = ("mean", "error", "variance", "tau", "bins", "jackknife")
    before = {field: np.array(getattr(value, field), copy=True) for field in fields}
    before_count = value.count
    with h5py.File(filename, "a") as archive:
        archive.copy("result", "broken")
        broken = archive["broken"]
        broken["count"][()] = 17
        broken["mean/value"][()] = before["mean"] + 100
        if corruption == "late_type":
            del broken["jacknife/data"]
            broken["jacknife/data"] = "not numeric"
        elif corruption == "missing_error":
            del broken["mean/error"]
        elif corruption == "zero_binsize":
            broken["timeseries/data"].attrs["binsize"] = np.uint64(0)
        elif corruption == "jack_count":
            del broken["jacknife/data"]
            broken["jacknife/data"] = before["jackknife"][:-1]
        elif corruption == "malformed_group":
            del broken["tau/value"]
            broken["tau/value/0"] = "not numeric"
        else:
            del broken["tau/value"]
            broken["tau/value"] = np.ones(3)
    with pytest.raises((RuntimeError, hdf5.ArchiveError)):
        value.load(filename, "/broken")
    assert value.count == before_count
    for field in fields:
        np.testing.assert_array_equal(getattr(value, field), before[field])
    # Metadata that the Python class does not expose must survive as well.
    value.save(filename, "/after_failure")
    with h5py.File(filename, "r") as archive:
        for field in ("binsize", "maxbinnum"):
            assert (archive["after_failure/timeseries/data"].attrs[field]
                    == archive["result/timeseries/data"].attrs[field])
