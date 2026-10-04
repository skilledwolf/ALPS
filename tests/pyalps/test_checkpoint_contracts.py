"""User workflows discovered by auditing the migration beyond API presence."""

import copy
import gc
import os
import subprocess
import sys
import textwrap
import weakref

import numpy as np
import pytest

from pyalps import hdf5, ngs


def result(vector=False, offset=0):
    observable = (ngs.createRealVectorObservable if vector else ngs.createRealObservable)("samples")
    for i in range(64):
        observable << (np.array([i + offset, 2.0 * i]) if vector else float(i + offset))
    return ngs.observable2result(observable)


def _measurement_samples(count, vector):
    indices = np.arange(count)
    scalar = (indices % 17 - 8) + 0.125 * (indices % 3)
    return np.column_stack((scalar, 2.0 * scalar + 3.0)) if vector else scalar


@pytest.mark.parametrize("vector", [False, True])
def test_measurement_checkpoint_continues_partial_bins(tmp_path, vector):
    samples = _measurement_samples(1031, vector)
    measurements = ngs.observables()
    create = (measurements.createRealVectorObservable if vector
              else measurements.createRealObservable)
    name = "energy"
    create(name)
    # An odd checkpoint boundary retains a partially filled bin and the
    # logarithmic binning state needed by subsequent samples.
    for sample in samples[:517]:
        measurements[name] << sample
    filename = str(tmp_path / "measurements.h5")
    with hdf5.archive(filename, "w") as archive:
        archive["custom/measurements"] = measurements
    restored = ngs.observables()
    with hdf5.archive(filename, "r") as archive:
        archive.set_context("/custom")
        restored.load(archive, "measurements")
        assert archive.context == "/custom"
    assert set(restored) == {name}
    assert ngs.observable2result(restored[name]).count == 517
    for sample in samples[517:]:
        measurements[name] << sample
        restored[name] << sample
    uninterrupted = ngs.observable2result(measurements[name])
    continued = ngs.observable2result(restored[name])
    assert continued.count == uninterrupted.count == len(samples)
    np.testing.assert_allclose(continued.mean, samples.mean(axis=0), rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(continued.variance, samples.var(axis=0, ddof=1),
                               rtol=1e-12, atol=1e-12)
    for attribute in ("mean", "error", "variance", "tau"):
        np.testing.assert_array_equal(getattr(continued, attribute),
                                      getattr(uninterrupted, attribute))
    # Vector result bins are available in the archive, although the Python
    # result.bins property currently supports only scalar results.
    with hdf5.archive(str(tmp_path / "continued-results.h5"), "w") as archive:
        archive["continued"] = continued
        archive["uninterrupted"] = uninterrupted
        np.testing.assert_array_equal(archive["continued/timeseries/data"],
                                      archive["uninterrupted/timeseries/data"])


@pytest.mark.parametrize("vector", [False, True])
def test_aligned_measurement_merge_error_and_result_checkpoint(tmp_path, vector):
    samples = _measurement_samples(256, vector)
    create = ngs.createRealVectorObservable if vector else ngs.createRealObservable
    left, right = create("energy"), create("energy")
    for sample in samples[:128]:
        left << sample
    for sample in samples[128:]:
        right << sample
    right_before = ngs.observable2result(right)
    left.merge(right)
    merged = ngs.observable2result(left)
    assert merged.count == len(samples)
    np.testing.assert_allclose(merged.mean, samples.mean(axis=0), rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(merged.variance, samples.var(axis=0, ddof=1),
                               rtol=1e-12, atol=1e-12)
    # Both runs end on complete bins. The merged result retains 128 bins,
    # each averaging two raw samples. Calculate its standard error directly
    # from the raw bin values, rather than another ALPS result.
    bins = samples.reshape((128, 2) + samples.shape[1:]).mean(axis=1)
    np.testing.assert_allclose(merged.error,
                               bins.std(axis=0, ddof=1) / np.sqrt(len(bins)),
                               rtol=1e-12, atol=1e-12)
    right_after = ngs.observable2result(right)
    assert right_after.count == right_before.count == 128
    for attribute in ("mean", "error", "variance", "tau"):
        np.testing.assert_array_equal(getattr(right_after, attribute),
                                      getattr(right_before, attribute))
    filename = str(tmp_path / "merged-results.h5")
    with hdf5.archive(filename, "w") as archive:
        archive["results/energy"] = merged
        archive["right_before"] = right_before
        archive["right_after"] = right_after
        np.testing.assert_array_equal(archive["right_before/timeseries/data"],
                                      archive["right_after/timeseries/data"])
    restored = ngs.result()
    with hdf5.archive(filename, "r") as archive:
        archive.set_context("/results/energy")
        restored.load(archive)
        assert archive.context == "/results/energy"
        np.testing.assert_allclose(archive["timeseries/data"], bins,
                                   rtol=1e-12, atol=1e-12)
    assert restored.count == merged.count
    assert np.asarray(restored.mean).shape == samples.shape[1:]
    for attribute in ("mean", "error", "variance", "tau"):
        np.testing.assert_array_equal(getattr(restored, attribute),
                                      getattr(merged, attribute))
    with hdf5.archive(filename, "a") as archive:
        archive["restored"] = restored
        np.testing.assert_array_equal(archive["restored/timeseries/data"],
                                      archive["results/energy/timeseries/data"])


def test_parameters_own_values_and_require_reassignment(tmp_path):
    array = np.array([1.0, 2.0])
    reference = weakref.ref(array)
    parameters = ngs.params({"array": array, "integers": [1, 2]})
    array[0] = 9
    del array
    gc.collect()
    assert reference() is None
    np.testing.assert_array_equal(parameters["array"], [1, 2])
    detached = parameters["array"]
    detached[0] = 10
    np.testing.assert_array_equal(parameters["array"], [1, 2])
    parameters["array"] = detached
    duplicate = copy.deepcopy(parameters)
    parameters["array"] += 1
    np.testing.assert_array_equal(parameters["array"], [11, 3])
    np.testing.assert_array_equal(duplicate["array"], [10, 2])


def test_parameter_checkpoint_retains_native_types_and_large_integers(tmp_path):
    values = {"array": np.arange(6.0), "large": 2 ** 53 + 1,
              "complex": 2 + 3j, "flags": [True, False], "empty": np.array([], dtype=bool)}
    filename = str(tmp_path / "parameters.h5")
    with hdf5.archive(filename, "w") as archive:
        archive["custom"] = ngs.params(values)
    with hdf5.archive(filename, "r") as archive:
        parameters = ngs.params(archive, "/custom")
        assert archive.context == "/"
    np.testing.assert_array_equal(parameters["array"], values["array"])
    assert parameters["large"] == 2 ** 53 + 1
    assert parameters["complex"] == 2 + 3j
    np.testing.assert_array_equal(parameters["flags"], [True, False])
    assert parameters["empty"].dtype.kind == "b"
    for unsupported in (np.arange(6).reshape(2, 3), {"label": "run"}, [True, 1], None):
        with pytest.raises(TypeError):
            parameters["unsupported"] = unsupported
    assert "unsupported" not in parameters


@pytest.mark.parametrize("vector", [False, True])
def test_result_loads_default_and_existing_objects_without_writing(tmp_path, vector):
    original = result(vector)
    filename = str(tmp_path / "result.h5")
    with hdf5.archive(filename, "w") as archive:
        archive["result"] = original
    restored = ngs.result()
    assert restored.count == 0 and repr(restored) == "No Measurements"
    with hdf5.archive(filename, "r") as archive:
        archive.set_context("/result")
        restored.load(archive)
        assert archive.context == "/result"
    np.testing.assert_array_equal(restored.mean, original.mean)
    np.testing.assert_array_equal(restored.error, original.error)
    assert restored.count == original.count
    alias = ngs.result(restored)
    replacement = result(not vector, offset=100)
    with hdf5.archive(filename, "w") as archive:
        archive["result"] = replacement
    with hdf5.archive(filename, "r") as archive:
        archive.set_context("/result")
        restored.load(archive)
    np.testing.assert_array_equal(alias.mean, original.mean)
    np.testing.assert_array_equal(restored.mean, replacement.mean)


def test_result_collection_load_uses_context_and_keeps_old_references(tmp_path):
    first = ngs.results()
    first["energy/with&name"] = result()
    first["vector"] = result(True)
    filename = str(tmp_path / "results.h5")
    with hdf5.archive(filename, "w") as archive:
        archive["custom/results"] = first
    restored = ngs.results()
    restored["old"] = result(offset=10)
    old_reference = restored["old"]
    with hdf5.archive(filename, "r") as archive:
        restored.load(archive, "/custom/results")
        assert archive.context == "/"
        assert set(restored) == set(first)
        energy_reference = restored["energy/with&name"]
        restored.load(archive, "/custom/results")
    del restored
    gc.collect()
    assert old_reference.mean == 41.5
    assert energy_reference.mean == 31.5


def test_failed_loads_preserve_values_and_context(tmp_path):
    filename = str(tmp_path / "broken.h5")
    with hdf5.archive(filename, "w") as archive:
        archive["bad/mean/value"] = 100.0  # no mandatory count dataset
    value = result()
    with hdf5.archive(filename, "r") as archive:
        archive.set_context("/bad")
        with pytest.raises(Exception):
            value.load(archive)
        assert archive.context == "/bad"
        assert value.mean == 31.5
        archive.set_context("/")
        parameters = ngs.params({"kept": [1, 2]})
        with pytest.raises(Exception):
            parameters.load(archive, "/missing")
        assert archive.context == "/"
        np.testing.assert_array_equal(parameters["kept"], [1, 2])
        results = ngs.results()
        results["kept"] = value
        with pytest.raises(Exception):
            results.load(archive, "/bad")
        assert archive.context == "/" and results["kept"].mean == 31.5


@pytest.mark.parametrize("vector", [False, True])
def test_result_unary_operations_and_deepcopy_do_not_mutate_original(vector):
    original = result(vector)
    mean = np.array(original.mean)
    positive, negative, duplicate = +original, -original, copy.deepcopy(original)
    duplicate += 10
    np.testing.assert_array_equal(original.mean, mean)
    np.testing.assert_array_equal(positive.mean, mean)
    np.testing.assert_array_equal(negative.mean, -mean)
    np.testing.assert_array_equal(duplicate.mean, mean + 10)


def test_empty_result_operations_raise_instead_of_crashing():
    code = '''
        from pyalps import ngs
        result = ngs.result()
        for operation in (lambda: result + 1, lambda: abs(result),
                          lambda: result ** 2, lambda: result.sin(),
                          lambda: result.mean, lambda: result.error):
            try:
                operation()
            except (RuntimeError, ValueError):
                pass
            else:
                raise AssertionError("empty result operation should fail")
    '''
    completed = subprocess.run(
        [sys.executable, "-X", "faulthandler", "-c", textwrap.dedent(code)],
        capture_output=True, text=True, timeout=30,
        env={**os.environ, "MallocScribble": "1"},
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr


def test_parameter_helpers_and_analysis_read_typed_checkpoints(tmp_path):
    import pyalps
    values = {"L": 8, "T": 1.5, "wide": 2 ** 53 + 1, "label": "a,b"}
    paths = pyalps.writeInputH5Files(str(tmp_path / "run"), [values])
    assert pyalps.getParameters(paths) == [values]
    props = pyalps.loadProperties(paths)[0]
    assert props["wide"] == 2 ** 53 + 1 and props["label"] == "a,b"
    with hdf5.archive(paths[0], "r") as archive:
        assert archive["/parameters/format"] == "alps.params.v2"


@pytest.mark.parametrize("vector", [False, True])
def test_analysis_loader_transfers_to_native_alea_for_results(tmp_path, vector):
    import pyalps

    value = result(vector)
    filename = str(tmp_path / "analysis.h5")
    with hdf5.archive(filename, "w") as archive:
        archive["parameters"] = ngs.params({"L": 2})
        archive["simulation/results/energy"] = value
    loaded = pyalps.loadMeasurements([filename], what=["energy"])
    assert len(loaded) == 1 and len(loaded[0]) == 1
    measured = loaded[0][0].y
    if vector:
        np.testing.assert_array_equal(measured.mean, value.mean)
        np.testing.assert_array_equal(measured.error, value.error)
    else:
        assert measured[0].mean == value.mean
        assert measured[0].error == value.error


def test_save_results_empty_collection_replaces_old_results(tmp_path):
    filename = str(tmp_path / "empty-results.h5")
    populated = ngs.results()
    populated["old"] = result()
    with hdf5.archive(filename, "w") as archive:
        ngs.saveResults(populated, ngs.params({"seed": 42}), archive, "/simulation/results")
        archive["unrelated"] = np.int64(9)
        ngs.saveResults(ngs.results(), ngs.params({"seed": 91}), archive, "/simulation/results")
        assert archive.is_group("/simulation/results")
        assert archive.list_children("/simulation/results") == []
        assert ngs.params(archive, "/parameters")["seed"] == 91
        assert archive["unrelated"] == 9
    restored = ngs.results()
    with hdf5.archive(filename, "r") as archive:
        restored.load(archive, "/simulation/results")
    assert len(restored) == 0


def test_extended_parameter_scalars_reject_lossy_conversion():
    if np.finfo(np.longdouble).nmant <= np.finfo(np.float64).nmant:
        pytest.skip("platform has no extended floating-point precision")
    value = np.longdouble(1) + np.finfo(np.longdouble).eps
    for supplied in (value, np.clongdouble(value + 1j)):
        with pytest.raises(TypeError, match="lossy"):
            ngs.params({"x": supplied})


def test_parameter_checkpoint_uses_native_hdf5_types_and_zero_extents(tmp_path):
    import h5py

    values = {"flag": True, "complex": 2 - 3j, "complexes": np.array([1 + 2j, 3 - 4j]),
              "empty flags": np.array([], dtype=bool),
              "empty complex": np.array([], dtype=np.complex128), "empty string": ""}
    filename = tmp_path / "native-params.h5"
    with hdf5.archive(str(filename), "w") as archive:
        archive["parameters"] = ngs.params(values)
    with h5py.File(filename, "r") as archive:
        assert archive["parameters/format"].asstr()[()] == "alps.params.v2"
        entries = {entry["name"].asstr()[()]: entry
                   for entry in archive["parameters/entries"].values()}
        assert entries["flag"]["value"].dtype == np.dtype(bool)
        assert entries["complex"]["value"].dtype == np.dtype(np.complex128)
        assert entries["complex"]["value"].shape == ()
        assert entries["complexes"]["value"].shape == (2,)
        assert entries["empty flags"]["value"].shape == (0,)
        assert entries["empty complex"]["value"].shape == (0,)
        assert entries["empty string"]["value"].asstr()[()] == ""
        for entry in entries.values():
            assert set(entry) == {"name", "value"}
            assert not any(name in entry["value"].attrs
                           for name in ("__complex__", "__alps_type__"))


@pytest.mark.parametrize("fault", ["v1", "null-vector", "pair-real-complex", "byte-bool", "extra-axis"])
def test_parameter_checkpoint_rejects_obsolete_or_mismatched_payloads(tmp_path, fault):
    import h5py

    filename = tmp_path / "bad-params.h5"
    original = ngs.params({"x": np.array([1 + 2j])})
    with hdf5.archive(str(filename), "w") as archive:
        archive["parameters"] = original
    with h5py.File(filename, "r+") as archive:
        group = archive["parameters"]
        entry = group["entries/0"]
        if fault == "v1":
            del group["format"]
            group["format"] = "alps.params.v1"
        else:
            del entry["value"]
            if fault == "null-vector":
                entry["value"] = h5py.Empty(np.complex128)
            elif fault == "pair-real-complex":
                value = entry.create_dataset("value", data=np.array([[1., 2.]]))
                value.attrs["__complex__"] = np.int8(1)
            elif fault == "byte-bool":
                entry["value"] = np.array([0, 1], dtype="i1")
            else:
                entry["value"] = np.array([[1 + 2j]])
    parameters = ngs.params({"kept": 42})
    with hdf5.archive(str(filename), "r") as archive:
        with pytest.raises(Exception):
            parameters.load(archive, "/parameters")
        assert archive.context == "/"
    assert parameters["kept"] == 42 and list(parameters) == ["kept"]
