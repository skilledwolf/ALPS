"""Check dtype-dependent user operations, not just equal numeric values."""

import h5py
import numpy as np
import pytest

from pyalps import hdf5, ngs


@pytest.mark.parametrize("depth", [2, 3])
@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize(
    "second_row,dtype",
    [
        ([np.int32(3), np.int32(4)], np.int32),
        ([np.float64(3), np.float64(4)], np.float64),
        ([np.complex128(3), np.complex128(4)], np.complex128),
        (np.array([3, 4], dtype=np.int32), np.int32),
    ],
)
def test_rectangular_python_numpy_rows_keep_array_contract(
    tmp_path, second_row, dtype, reverse, depth
):
    # Both rows have the same ALPS storage dtype. Plain Python integers use
    # int32, even when NumPy's platform-default integer is int64.
    scalar = {np.int32: int, np.float64: float, np.complex128: complex}[dtype]
    value = [[scalar(1), scalar(2)], second_row]
    if reverse:
        value.reverse()
    if depth == 3:
        value = [value, value]
    expected = np.asarray(value, dtype=dtype)
    filename = str(tmp_path / "rectangular.h5")
    with hdf5.archive(filename, "w") as archive:
        archive["table"] = value
    with hdf5.archive(filename, "r") as archive:
        restored = archive["table"]
        assert archive.is_data("table")
    # Value-only NumPy comparisons coerce lists and miss changed arithmetic.
    assert isinstance(restored, np.ndarray)
    assert restored.dtype == expected.dtype and restored.shape == expected.shape
    np.testing.assert_array_equal(restored[..., 0], expected[..., 0])
    np.testing.assert_array_equal(restored * 2, expected * 2)


@pytest.mark.parametrize(
    "value",
    [
        [[1, np.int32(2)], [np.int32(3), np.int32(4)]],
        [[1, 2], [np.int32(3)]],
        [[True, False], [np.int32(3), np.int32(4)]],
        [[2**53 + 1, 2**53 + 3], [np.int32(3), np.int32(4)]],
    ],
)
def test_incompatible_python_numpy_rows_keep_groups(tmp_path, value):
    with hdf5.archive(str(tmp_path / "groups.h5"), "w") as archive:
        archive["table"] = value
        assert archive.is_group("table")
        restored = archive["table"]
    assert isinstance(restored, list)
    for actual, expected in zip(restored, value):
        assert list(actual) == list(expected)


@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("nested", [False, True])
@pytest.mark.parametrize(
    "other_dtype,integer_dtype,integers",
    [
        (np.int64, np.uint64, [2**63 + 1, 2**63 + 3]),
        (np.uint64, np.int64, [2**53 + 1, 2**53 + 3]),
        (np.float32, np.int64, [2**53 + 1, 2**53 + 3]),
        (np.float64, np.int64, [-(2**53) - 1, -(2**53) - 3]),
        (np.complex64, np.uint64, [2**64 - 1, 2**64 - 3]),
        (np.complex128, np.int64, [2**63 - 1, -(2**63) + 1]),
        (None, np.int64, [2**53 + 1, 2**53 + 3]),
    ],
)
def test_lossy_mixed_rows_preserve_exact_integers(
    tmp_path, other_dtype, integer_dtype, integers, reverse, nested
):
    other = np.array([3, 4], dtype=other_dtype) if other_dtype else [3.0, 4.0]
    # Exercise a noncontiguous, read-only integer row too.
    row = np.repeat(np.array(integers, dtype=integer_dtype), 2)[::2]
    row.flags.writeable = False
    rows = [other, row]
    if reverse:
        rows.reverse()
    value = [rows, rows] if nested else rows
    filename = str(tmp_path / "exact-integers.h5")
    with hdf5.archive(filename, "w") as archive:
        archive["table"] = np.zeros((2, 2))  # Replace an existing dataset.
        archive["table"] = value
    with hdf5.archive(filename, "r") as archive:
        assert archive.is_group("table")
        restored = archive["table"]
    for actual in restored if nested else [restored]:
        integer_row = actual[0 if reverse else 1]
        assert integer_row.dtype == integer_dtype
        # NumPy mixed int/float equality would silently round the expected
        # integers too. Compare independently converted Python ints.
        assert [int(x) for x in integer_row] == integers
        np.testing.assert_array_equal(actual[1 if reverse else 0], other)


@pytest.mark.parametrize("scalar", [int, np.int64])
@pytest.mark.parametrize("nested", [False, True])
def test_integer_scalar_rows_are_checked_before_stacking(tmp_path, scalar, nested):
    integers = [2**53 + 1, 2**53 + 3]
    rows = [np.array([1.5, 2.5]), tuple(scalar(x) for x in integers)]
    value = [rows] if nested else rows
    filename = str(tmp_path / "scalar-rows.h5")
    with hdf5.archive(filename, "w") as archive:
        archive["table"] = value
    with hdf5.archive(filename, "r") as archive:
        restored = archive["table"]
    if nested:
        restored = restored[0]
    assert [int(x) for x in restored[1]] == integers


@pytest.mark.parametrize(
    "other_dtype", [np.uint64, np.float32, np.float64, np.complex64]
)
@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("integers", [[0, 2**40], [2**53, 2**53 + 2], [-(2**63), 0]])
def test_lossless_mixed_rows_still_stack(tmp_path, other_dtype, reverse, integers):
    other = np.array([3, 4], dtype=other_dtype)
    rows = [other, np.array(integers, dtype=np.int64)]
    if reverse:
        rows.reverse()
    filename = str(tmp_path / "lossless.h5")
    with hdf5.archive(filename, "w") as archive:
        archive["table"] = {"old": 1}  # Replace an existing group.
        archive["table"] = rows
    with hdf5.archive(filename, "r") as archive:
        assert archive.is_data("table")
        restored = archive["table"]
    assert isinstance(restored, np.ndarray)
    assert restored.shape == (2, 2)
    assert restored.dtype == np.asarray(rows).dtype
    assert [int(x.real) for x in restored[0 if reverse else 1]] == integers


def test_lossless_mixed_rows_allow_nan_and_infinity(tmp_path):
    with hdf5.archive(str(tmp_path / "nonfinite.h5"), "w") as archive:
        archive["table"] = [np.array([np.nan, np.inf]), np.array([2**53, 0])]
        assert archive.is_data("table")
        restored = archive["table"]
    assert np.isnan(restored[0, 0]) and np.isposinf(restored[0, 1])
    assert [int(x) for x in restored[1]] == [2**53, 0]


@pytest.mark.parametrize("shape", [(), (0,), (2, 0)])
def test_mixed_integer_array_shapes(tmp_path, shape):
    rows = [np.zeros(shape, dtype=np.int64), np.full(shape, 2**64 - 1, dtype=np.uint64)]
    filename = str(tmp_path / "shapes.h5")
    with hdf5.archive(filename, "w") as archive:
        archive["table"] = rows
    with hdf5.archive(filename, "r") as archive:
        restored = archive["table"]
    if shape == ():
        assert restored == [0, 2**64 - 1]
    else:
        assert isinstance(restored, np.ndarray)
        assert restored.shape == (2, *shape)


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
            assert type(actual) is (bool if dtype == np.bool_ else int)
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
        assert archive[path] is value


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
        assert type(actual) is {np.bool_: bool, np.int8: int, np.complex128: complex}[dtype]
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
    native = value.astype(value.dtype.newbyteorder("="))
    assert actual.dtype == native.dtype
    assert actual.real.tobytes() == native.real.tobytes()
    assert actual.imag.tobytes() == native.imag.tobytes()


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
        assert archive.extent(path) == []
        with pytest.raises(hdf5.WrongType, match="NULL.*shape.*offline converter"):
            archive[path]
