import gc
import tempfile
import weakref
import numpy as np
from pyalps import hdf5, ngs
import parameter_probe as native

empty = native.empty_vectors()
for key, kind in {"integer": "i", "real": "f", "complex": "c", "boolean": "b"}.items():
    assert empty[key].size == 0 and empty[key].dtype.kind == kind
assert empty["text"] == []

parameters = native.string_parameters()
names = parameters["value"]
names[0] = "Magnetization"
assert native.string_vector(parameters) == ["Sz"]
parameters["value"] = names
assert native.string_vector(parameters) == ["Magnetization"]

values = np.array([1., 2.])
reference = weakref.ref(values)
parameters = ngs.params({"value": values})
values *= 3
assert native.vector(parameters) == [1., 2.]
del values
gc.collect()
assert reference() is None
copy = parameters["value"]
copy[0] = 8
assert native.threaded_vector(parameters) == [1., 2.]
parameters["value"] = copy
clone = native.clone(parameters)
native.replace(parameters)
assert native.vector(parameters) == [5., 6.]
assert native.vector(clone) == [8., 2.]
parameters["value"] = 2 ** 53 + 1
assert native.wide_integer(parameters) == 2 ** 53 + 1
for reader in (native.integer, native.vector):
    try:
        reader(parameters)
    except (RuntimeError, TypeError, ValueError, IndexError, OverflowError):
        pass
    else:
        raise AssertionError("incompatible or narrowing conversion accepted")
parameters["value"] = [1, 2.5]
assert native.vector(parameters) == [1., 2.5]
assert native.complex_vector(parameters) == [1+0j, 2.5+0j]
try:
    native.integer_vector(parameters)
except (RuntimeError, TypeError, ValueError):
    pass
else:
    raise AssertionError("real to integer conversion accepted")

with tempfile.TemporaryDirectory() as directory:
    filename = directory + "/shared-archive.h5"
    with hdf5.archive(filename, "w") as archive:
        archive["initial"] = 1
    # Each provider owns its own file. Close h5py before a native writer opens it.
    native.append_archive(filename)
    with hdf5.archive(filename, "r") as archive:
        assert archive["initial"] == 1 and archive["from_native"] == 42
    for i, values in enumerate(([True, False], [2**53+1, 2], [1+2j, 3+4j], ["a,b", "c"])):
        with hdf5.archive(directory + f"/params-{i}.h5", "w") as archive:
            archive["parameters"] = ngs.params({"value": values})
            archive.set_context("/parameters")
            restored = native.empty_vectors()
            with archive.native() as native_archive:
                native.load(restored, native_archive)
            np.testing.assert_array_equal(restored["value"], values)
            assert archive.context == "/parameters"
            archive.set_context("/resaved")
            with archive.native() as native_archive:
                native.save(restored, native_archive)
                native.load(restored, native_archive)
            np.testing.assert_array_equal(restored["value"], values)
            assert archive.context == "/resaved"
    # Legacy parameter groups are explicitly rejected, without partial load.
    with hdf5.archive(directory + "/legacy.h5", "w") as archive:
        archive["value"] = [2., 4.]
        try:
            with archive.native() as native_archive:
                native.load(clone, native_archive)
        except RuntimeError:
            pass
        else:
            raise AssertionError("legacy checkpoint accepted")
        assert native.threaded_vector(clone) == [8., 2.]
native.destroy_on_worker(parameters)
print("native parameter contracts: ok")
