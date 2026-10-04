import numpy as np
import pyalps.hdf5 as hdf5


# 2**20 doubles = 8 MiB, enough to exercise array IO without making the
# suite slow. The original's 2**29 ceiling tested the operating system, not ALPS.
MAX_EXPONENT = 20


def test_growing_arrays_roundtrip(tmp_path):
    """Explicit array writes retain their dtype and values."""
    pattern = str(tmp_path / "foo.h5")

    written = {}
    with hdf5.archive(pattern, "a") as archive:
        size = 2 ** 10
        while size <= 2 ** MAX_EXPONENT:
            values = np.linspace(0.0, 1.0, size)
            archive[str(size)] = values
            written[str(size)] = values
            size *= 2

    with hdf5.archive(pattern, "r") as archive:
        for key, values in written.items():
            restored = archive[key]
            assert restored.shape == values.shape, key
            assert restored.dtype == values.dtype, key
            np.testing.assert_array_equal(restored, values)


def test_empty_and_single_element_arrays(tmp_path):
    """The size-0 and size-1 edges of the same write path."""
    path = str(tmp_path / "edges.h5")

    with hdf5.archive(path, "a") as archive:
        archive["/empty"] = np.empty(0)
        archive["/single"] = np.array([np.pi])

    with hdf5.archive(path, "r") as archive:
        assert archive["/empty"].shape == (0,)
        np.testing.assert_array_equal(archive["/single"], np.array([np.pi]))
