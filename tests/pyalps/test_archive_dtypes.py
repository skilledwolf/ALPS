"""Check dtype-dependent user operations, not just equal numeric values."""

import h5py
import numpy as np
import pytest

from pyalps import hdf5, ngs


@pytest.mark.parametrize("dtype", [np.bool_, np.int8])
@pytest.mark.parametrize("shape", [(), (3,), (2, 3), (2, 1, 3), (0,), (2, 0)])
@pytest.mark.parametrize("attribute", [False, True])
def test_boolean_and_signed_byte_round_trip(tmp_path, dtype, shape, attribute):
    size = int(np.prod(shape))
    value = np.arange(size, dtype=np.int8).reshape(shape).astype(dtype)
    if value.ndim > 1:
        value = np.asfortranarray(value)
    value.flags.writeable = False
    path = "/group/@value" if attribute else "/value"
    filename = str(tmp_path / "dtype.h5")
    with hdf5.archive(filename, "w") as archive:
        archive.create_group("/group")
        archive[path] = value
    with hdf5.archive(filename, "r") as archive:
        actual = archive[path]
        if shape:
            assert actual.dtype == dtype
            assert actual.shape == shape
        else:
            assert np.asarray(actual).dtype == dtype
        np.testing.assert_array_equal(actual, value)


@pytest.mark.parametrize("scalar", [bool, np.bool_])
@pytest.mark.parametrize("value", [False, True])
@pytest.mark.parametrize("attribute", [False, True])
def test_boolean_scalar_round_trip(tmp_path, scalar, value, attribute):
    path = "/group/@flag" if attribute else "/flag"
    filename = str(tmp_path / "boolean-scalar.h5")
    with hdf5.archive(filename, "w") as archive:
        archive.create_group("/group")
        archive[path] = scalar(value)
    with hdf5.archive(filename, "r") as archive:
        assert bool(archive[path]) is value


@pytest.mark.parametrize("value", [False, True])
def test_numpy_boolean_scalar_parameter_checkpoint(tmp_path, value):
    filename = str(tmp_path / "boolean-parameter.h5")
    parameters = ngs.params({"flag": np.bool_(value)})
    with hdf5.archive(filename, "w") as archive:
        archive["parameters"] = parameters
    with hdf5.archive(filename, "r") as archive:
        restored = ngs.params(archive, "/parameters")
    assert restored["flag"] is value


def test_boolean_mask_remains_a_mask_after_reload(tmp_path):
    with hdf5.archive(str(tmp_path / "mask.h5"), "w") as archive:
        archive["mask"] = np.array([True, False, True])
        np.testing.assert_array_equal(np.array([10, 20, 30])[archive["mask"]], [10, 30])
        # Cover the native vector<bool> writer as well as ndarray dispatch.
        archive["parameters"] = ngs.params({"mask": [True, False, True]})
        np.testing.assert_array_equal(
            np.array([10, 20, 30])[ngs.params(archive, "/parameters")["mask"]], [10, 30]
        )


@pytest.mark.parametrize("attribute", [False, True])
def test_overwriting_boolean_and_integer_replaces_storage_type(tmp_path, attribute):
    path = "/group/@value" if attribute else "/value"
    with hdf5.archive(str(tmp_path / "overwrite.h5"), "w") as archive:
        archive.create_group("/group")
        for dtype in (np.bool_, np.int8, np.bool_, np.int8):
            archive[path] = np.array([0, 1], dtype=dtype)
            assert archive[path].dtype == dtype


@pytest.mark.parametrize("attribute", [False, True])
def test_unmarked_signed_byte_is_integer_even_for_zero_one_values(tmp_path, attribute):
    filename = str(tmp_path / "signed-byte.h5")
    path = "/@value" if attribute else "/value"
    with h5py.File(filename, "w") as archive:
        if attribute:
            archive.attrs["value"] = np.array([0, 1, 0], dtype=np.int8)
        else:
            archive["value"] = np.array([0, 1, 0], dtype=np.int8)
    with hdf5.archive(filename, "r") as archive:
        actual = archive[path]
    assert actual.dtype == np.int8
    np.testing.assert_array_equal(actual + 2, [2, 3, 2])


@pytest.mark.parametrize("dtype", [np.bool_, np.int8, np.complex64, np.complex128])
@pytest.mark.parametrize("shape", [(), (3,), (2, 3), (0,), (2, 0)])
@pytest.mark.parametrize("attribute", [False, True])
def test_native_storage_uses_h5py_datatypes_without_private_markers(
    tmp_path, dtype, shape, attribute
):
    size = int(np.prod(shape))
    value = np.arange(size).reshape(shape).astype(dtype)
    if np.issubdtype(dtype, np.complexfloating):
        value.imag = -np.arange(size).reshape(shape)
    filename = str(tmp_path / "standard-types.h5")
    path = "/group/@value" if attribute else "/value"
    with hdf5.archive(filename, "w") as archive:
        archive.create_group("/group")
        archive[path] = value
        assert archive.extent(path) == list(shape)
        assert archive.is_scalar(path) == (shape == ())
        assert archive.is_complex(path) == np.issubdtype(dtype, np.complexfloating)
    with h5py.File(filename, "r") as archive:
        parent = archive["group"] if attribute else archive["value"]
        actual = parent.attrs["value"] if attribute else parent[()]
        datatype = parent.attrs.get_id("value").get_type() if attribute else parent.id.get_type()
        try:
            assert datatype.get_class() == (
                h5py.h5t.ENUM if dtype == np.bool_
                else h5py.h5t.COMPOUND if np.issubdtype(dtype, np.complexfloating)
                else h5py.h5t.INTEGER
            )
            if dtype == np.bool_:
                assert {datatype.get_member_name(i): datatype.get_member_value(i)
                        for i in range(datatype.get_nmembers())} == {b"FALSE": 0, b"TRUE": 1}
            elif np.issubdtype(dtype, np.complexfloating):
                assert [datatype.get_member_name(i) for i in range(datatype.get_nmembers())] == [b"r", b"i"]
        finally:
            datatype.close()
        assert np.asarray(actual).dtype == dtype
        assert np.asarray(actual).shape == shape
        np.testing.assert_array_equal(actual, value)
        assert not any(name.startswith(("__complex__", "__alps_type__")) for name in parent.attrs)


@pytest.mark.parametrize("dtype", [np.bool_, np.int8, np.complex64, np.complex128])
@pytest.mark.parametrize("shape", [(), (3,), (2, 3), (0,), (2, 0)])
@pytest.mark.parametrize("attribute", [False, True])
def test_h5py_standard_types_load_without_alps_metadata(tmp_path, dtype, shape, attribute):
    size = int(np.prod(shape))
    value = np.arange(size).reshape(shape).astype(dtype)
    if np.issubdtype(dtype, np.complexfloating):
        value.imag = -np.arange(size).reshape(shape)
    filename = str(tmp_path / "from-h5py.h5")
    path = "/@value" if attribute else "/value"
    with h5py.File(filename, "w") as archive:
        if attribute:
            archive.attrs["value"] = value
        else:
            archive["value"] = value
    with hdf5.archive(filename, "r") as archive:
        actual = archive[path]
        assert archive.extent(path) == list(shape)
        assert archive.is_scalar(path) == (shape == ())
    if shape:
        assert actual.dtype == dtype
        assert actual.shape == shape
    elif dtype == np.complex64:
        assert np.asarray(actual).dtype == np.complex64
    else:
        assert np.asarray(actual).dtype == dtype
    np.testing.assert_array_equal(actual, value)


@pytest.mark.parametrize("dtype", ["<c8", ">c8", "<c16", ">c16"])
@pytest.mark.parametrize("attribute", [False, True])
def test_h5py_complex_preserves_endian_precision_and_ieee_components(tmp_path, dtype, attribute):
    value = np.empty(4, dtype=dtype)
    value.real = [0.0, np.inf, np.nan, -4.5]
    value.imag = [-0.0, -np.inf, 1.25, np.nan]
    filename = str(tmp_path / "ieee-complex.h5")
    path = "/@value" if attribute else "/value"
    with h5py.File(filename, "w") as archive:
        if attribute:
            archive.attrs["value"] = value
        else:
            archive["value"] = value
    with hdf5.archive(filename, "r") as archive:
        actual = archive[path]
    assert actual.dtype == value.dtype
    assert actual.real.tobytes() == value.real.tobytes()
    assert actual.imag.tobytes() == value.imag.tobytes()


@pytest.mark.parametrize("dtype", [np.bool_, np.int8, np.complex64, np.complex128])
@pytest.mark.parametrize("attribute", [False, True])
def test_null_dataspace_requires_explicit_shape_migration(tmp_path, dtype, attribute):
    filename = str(tmp_path / "null.h5")
    path = "/@value" if attribute else "/value"
    with h5py.File(filename, "w") as archive:
        if attribute:
            archive.attrs["value"] = h5py.Empty(dtype)
        else:
            archive.create_dataset("value", data=h5py.Empty(dtype))
    with hdf5.archive(filename, "r") as archive:
        assert archive.is_null(path)
        assert archive.extent(path) is None
        assert isinstance(archive[path], h5py.Empty)
        with archive.native() as native:
            with pytest.raises(hdf5.WrongType, match="NULL.*shape.*offline converter"):
                native[path]
