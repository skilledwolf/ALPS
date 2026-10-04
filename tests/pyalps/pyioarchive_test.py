"""Python files use ordinary h5py ownership and exact r/a/w modes."""
import pytest
from pyalps import hdf5


def test_archive_modes(tmp_path):
    filename = tmp_path / "modes.h5"
    with hdf5.archive(filename, "w") as archive:
        archive["old"] = 1
    with hdf5.archive(filename, "a") as archive:
        archive["new"] = 2
        assert archive["old"] == 1
    with hdf5.archive(filename, "w") as archive:
        assert archive.list_children("/") == []
        archive["fresh"] = 3
    with hdf5.archive(filename, "r") as archive:
        assert archive["fresh"] == 3
        with pytest.raises((OSError, RuntimeError, ValueError)):
            archive["fresh"] = 4
    for mode in ("al", "rl", "ac", "m", "", "rw"):
        with pytest.raises(ValueError, match="mode"):
            hdf5.archive(filename, mode)
