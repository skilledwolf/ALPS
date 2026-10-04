"""h5py primitive IO and native scientific checkpoint ownership boundaries."""
import h5py
import numpy as np
import pytest
from pyalps import hdf5, ngs


@pytest.mark.parametrize("dtype", ["i4", "i8", "u8", "f4", "f8", "c8", "c16", "bool"])
@pytest.mark.parametrize("shape", [(), (3,), (2, 3), (0,), (2, 0)])
def test_explicit_numpy_dtype_shape_and_noncontiguous_arrays(tmp_path, dtype, shape):
    value = np.arange(int(np.prod(shape))).reshape(shape).astype(dtype)
    if len(shape) == 2:
        value = np.asfortranarray(value)
    filename = tmp_path / "data.h5"
    with hdf5.archive(filename, "w") as archive:
        archive["value"] = value
        actual = archive["value"]
        assert np.asarray(actual).dtype == value.dtype
        assert np.asarray(actual).shape == value.shape
        np.testing.assert_array_equal(actual, value)
        assert archive.extent("value") == list(shape)
    with h5py.File(filename, "r") as independent:
        assert independent["value"].dtype == value.dtype
        np.testing.assert_array_equal(independent["value"][()], value)


def test_groups_have_no_inferred_python_container_type(tmp_path):
    with hdf5.archive(tmp_path / "groups.h5", "w") as archive:
        archive.create_group("items")
        archive["items/0"] = np.array([1, 2], dtype="i8")
        archive["items/1"] = np.array([3], dtype="i8")
        group = archive["items"]
        assert isinstance(group, h5py.Group)
        assert list(group) == ["0", "1"]
        with pytest.raises((ValueError, TypeError)):
            archive["items"] = [[1, 2], [3]]
        np.testing.assert_array_equal(archive["items/0"], [1, 2])


@pytest.mark.parametrize("dtype", ["i2", "i8", "u1", "u8", "f4", "f8", "c8", "c16", "bool"])
@pytest.mark.parametrize("shape", [(), (2, 3), (2, 0)])
def test_native_callbacks_preserve_numeric_scalar_dtype_and_array_rank(tmp_path, dtype, shape):
    value = np.arange(int(np.prod(shape))).reshape(shape).astype(dtype)
    with hdf5.archive(tmp_path / "native-dtype.h5", "w") as archive:
        with archive.native() as native:
            native["value"] = value
            actual = native["value"]
            assert np.asarray(actual).dtype == value.dtype
            assert np.asarray(actual).shape == shape
            if not shape:
                assert isinstance(actual, np.generic)
            np.testing.assert_array_equal(actual, value)
        np.testing.assert_array_equal(archive["value"], value)


@pytest.mark.parametrize("value", [np.array("é"),
                                  np.array([["", "é"], ["spin", "longer"]]),
                                  np.empty((2, 0), dtype="U8")])
def test_native_callbacks_preserve_unicode_array_shape(tmp_path, value):
    with hdf5.archive(tmp_path / "native-strings.h5", "w") as archive:
        with archive.native() as native:
            native["labels"] = value
            actual = native["labels"]
            assert np.asarray(actual).shape == value.shape
            assert np.asarray(actual).dtype.kind == "U"
            np.testing.assert_array_equal(actual, value)
        np.testing.assert_array_equal(archive["labels"], value)


def test_native_checkpoints_transfer_ownership_and_preserve_context(tmp_path):
    filename = tmp_path / "checkpoint.h5"
    parameters = ngs.params({"count": 2**53 + 1, "mask": np.array([True, False])})
    rng = ngs.random01(91)
    with hdf5.archive(filename, "w") as archive:
        archive["primitive"] = np.array([1.5, 2.5], dtype="f4")
        borrowed_group = archive["/"]
        archive.set_context("/state")
        archive["parameters"] = parameters
        assert not borrowed_group.id.valid  # h5py owner closed during native save
        assert archive.context == "/state"
        archive["random"] = rng
        restored = ngs.params(archive, "/state/parameters")
        assert restored["count"] == 2**53 + 1
        np.testing.assert_array_equal(restored["mask"], [True, False])
        with archive.native() as native:
            assert native.context == "/state"
            native["explicit-array"] = np.array([1, 2], dtype="i4")
            with pytest.raises(TypeError, match="explicit NumPy"):
                native["implicit-list"] = [1, 2]
            with pytest.raises(hdf5.ArchiveError, match="owns"):
                archive["primitive"]
            with pytest.raises(hdf5.ArchiveError, match="owns"):
                archive.set_context("/wrong-owner")
        np.testing.assert_array_equal(archive["explicit-array"], [1, 2])
        np.testing.assert_array_equal(archive["/primitive"], [1.5, 2.5])
        with pytest.raises((hdf5.ArchiveError, RuntimeError)):
            ngs.params(archive, "/missing")
        archive["after-error"] = np.int64(7)
        assert archive.context == "/state"
        assert archive["after-error"] == 7


def test_retained_native_callback_archive_is_closed_before_h5py_reopens(tmp_path):
    retained = []

    class Field:
        def save(self, native):
            native["value"] = np.array([3, 4], dtype="i4")

    class Simulation(ngs.mcbase):
        def update(self): pass
        def measure(self): pass
        def fraction_completed(self): return 1.0
        def save(self, native):
            retained.append(native)
            super().save(native)
            native["weights"] = np.array([1.5, 2.5], dtype="f8")
            native["field"] = Field()

    simulation = Simulation({"SEED": 42})
    simulation.measurements << ngs.RealObservable("energy")
    simulation.measurements["energy"] << 1.0
    with hdf5.archive(tmp_path / "retained.h5", "w") as archive:
        with archive.native() as native:
            native["simulation"] = simulation
            assert retained[0].is_open
        assert not retained[0].is_open
        with pytest.raises(hdf5.ArchiveClosed):
            retained[0]["weights"]
        np.testing.assert_array_equal(archive["simulation/weights"], [1.5, 2.5])
        np.testing.assert_array_equal(archive["simulation/field/value"], [3, 4])
        archive["after-native"] = np.int64(7)


def test_python_checkpoint_hook_uses_h5py_and_native_methods_serially(tmp_path):
    class State:
        def save(self, archive):
            archive["weights"] = np.array([1.25, 2.5], dtype="f8")
            archive["parameters"] = ngs.params({"seed": 91})
            archive["complete"] = np.bool_(True)

    with hdf5.archive(tmp_path / "state.h5", "w") as archive:
        archive["nested/state"] = State()
        assert archive.context == "/"
        np.testing.assert_array_equal(archive["nested/state/weights"], [1.25, 2.5])
        assert ngs.params(archive, "/nested/state/parameters")["seed"] == 91
        assert archive["nested/state/complete"]


def test_explicit_native_close_inside_transfer_reopens_h5py(tmp_path):
    with hdf5.archive(tmp_path / "closed-callback.h5", "w") as archive:
        with archive.native() as native:
            native["value"] = np.float32(1.5)
            native.close()
        assert archive["value"] == np.float32(1.5)
        archive["after-close"] = np.int64(7)
