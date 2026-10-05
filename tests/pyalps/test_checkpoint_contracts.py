"""User workflows discovered by auditing the migration beyond API presence."""

import copy
import gc
import weakref

import numpy as np
import pytest

from pyalps import alea, hdf5, ngs


def result(vector=False, offset=0):
    accumulator = alea.BatchAccumulator(2 if vector else 1)
    for i in range(64):
        accumulator << (np.array([i + offset, 2.0 * i]) if vector else float(i + offset))
    return accumulator.result()


class Simulation(ngs.mcbase):
    def __init__(self, parameters):
        super().__init__(parameters)
        self.measurements["energy/with&name"] = alea.BatchAccumulator(2, num_batches=8, base_size=3)
        self.steps = 0

    def update(self):
        self.sample = [self.random(), self.steps]
        self.steps += 1

    def measure(self):
        self.measurements["energy/with&name"] << self.sample

    def fraction_completed(self):
        return self.steps / 37

    def save(self, archive):
        super().save(archive)
        archive["checkpoint/steps"] = self.steps

    def load(self, archive):
        super().load(archive)
        self.steps = int(archive["checkpoint/steps"])


def test_mcbase_checkpoint_continues_rng_and_partial_native_batches(tmp_path):
    full, stopped = Simulation({"SEED": 7}), Simulation({"SEED": 7})
    full.run(lambda: False)
    for _ in range(11):
        stopped.update()
        stopped.measure()
    filename = tmp_path / "mcbase.h5"

    def save(archive):
        archive.set_context("/custom")
        stopped.save(archive)

    hdf5.save_checkpoint(filename, save)
    resumed = Simulation({"SEED": 99})
    old_handle = resumed.measurements["energy/with&name"]
    with hdf5.archive(filename) as archive:
        archive.set_context("/custom")
        resumed.load(archive)
        assert archive.context == "/custom"
        assert archive["measurements/energy&#47;with&#38;name/@kind"] == 6
    assert old_handle.count == 0
    assert resumed.run(lambda: False)
    expected, actual = ngs.collectResults(full), ngs.collectResults(resumed)
    for name in expected:
        assert actual[name].count == expected[name].count == 37
        for field in ("mean", "error", "covariance", "batch_sums", "batch_counts"):
            np.testing.assert_array_equal(getattr(actual[name], field), getattr(expected[name], field))
    assert [resumed.random() for _ in range(31)] == [full.random() for _ in range(31)]


@pytest.mark.parametrize("fault", ["cursor", "engine", "missing"])
def test_mcbase_failed_checkpoint_load_preserves_all_base_state(tmp_path, fault):
    import h5py
    saved = Simulation({"SEED": 7, "label": "new"})
    saved.update()
    saved.measure()
    filename = tmp_path / "broken-mcbase.h5"
    hdf5.save_checkpoint(filename, saved.save)
    with h5py.File(filename, "a") as archive:
        if fault == "engine":
            archive["checkpoint/engine/engine"][()] = "invalid random state"
        elif fault == "cursor":
            archive["measurements/energy&#47;with&#38;name/cursor/level"][()] = np.uint64(64)
        else:
            del archive["measurements/energy&#47;with&#38;name"]
    existing = Simulation({"SEED": 19, "label": "kept"})
    existing.update()
    existing.measure()
    handle = existing.measurements["energy/with&name"]
    before, random = handle.result(), copy.deepcopy(existing.random)
    with hdf5.archive(filename) as archive:
        with pytest.raises((RuntimeError, ValueError)):
            existing.load(archive)
        assert archive.context == "/"
    assert existing.parameters["label"] == "kept" and existing.steps == 1
    assert existing.measurements["energy/with&name"] is handle
    np.testing.assert_array_equal(handle.result().batch_sums, before.batch_sums)
    assert [existing.random() for _ in range(31)] == [random() for _ in range(31)]


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


def test_failed_parameter_load_preserves_values_and_context(tmp_path):
    filename = str(tmp_path / "broken.h5")
    with hdf5.archive(filename, "w") as archive:
        archive["bad/mean/value"] = 100.0
    parameters = ngs.params({"kept": [1, 2]})
    with hdf5.archive(filename, "r") as archive:
        archive.set_context("/bad")
        with pytest.raises(Exception):
            parameters.load(archive, "/missing")
        assert archive.context == "/bad"
        np.testing.assert_array_equal(parameters["kept"], [1, 2])


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
    np.testing.assert_array_equal([component.mean for component in measured], value.mean)
    np.testing.assert_array_equal([component.error for component in measured], value.error)


def test_save_results_empty_collection_replaces_old_results(tmp_path):
    filename = str(tmp_path / "empty-results.h5")
    populated = {}
    populated["old"] = result()
    with hdf5.archive(filename, "w") as archive:
        ngs.saveResults(populated, ngs.params({"seed": 42}), archive, "/simulation/results")
        archive["unrelated"] = np.int64(9)
        ngs.saveResults({}, ngs.params({"seed": 91}), archive, "/simulation/results")
        assert archive.is_group("/simulation/results")
        assert archive.list_children("/simulation/results") == []
        assert ngs.params(archive, "/parameters")["seed"] == 91
        assert archive["unrelated"] == 9


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
