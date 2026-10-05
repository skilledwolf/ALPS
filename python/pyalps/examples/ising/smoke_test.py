#!/usr/bin/env python3
"""Exercise the public downstream simulation-export helper with native ALEA batches."""

import faulthandler
import copy
import os
import shutil
import tempfile
import h5py
import numpy as np

# Report the Python frame if a native callback or import stalls in CI.
faulthandler.dump_traceback_later(30, exit=True)

# Initialize the package's library search directories before loading a
# downstream extension that shares its native runtime.
import pyalps
import ising_c
import pyalps.hdf5 as hdf5
import pyalps.ngs as ngs


parameters = ngs.params({"SEED": 7, "SWEEPS": 10,
                         "couplings": np.array([1., 2.]),
                         "label": "checkpoint"})
simulation = ising_c.sim(parameters)

assert issubclass(ising_c.sim, ngs.mcbase)
assert isinstance(simulation, ngs.mcbase)
assert int(simulation.parameters["SWEEPS"]) == 10
assert len(simulation.measurements) == 1
assert 0.0 <= simulation.random() < 1.0
# The base descriptors must work for a downstream C++ simulation too. This
# used to throw std::bad_cast because mcbase assumed every instance was its
# Python trampoline alias.
assert ngs.mcbase.parameters.__get__(simulation) is simulation.parameters
assert ngs.mcbase.measurements.__get__(simulation) is simulation.measurements
assert ngs.mcbase.random.__get__(simulation) is simulation.random
assert simulation.run(lambda: False)
assert simulation.resultNames() == ["Magnetization"]
before = simulation.collectResults()
assert before["Magnetization"].count == 10
simulation.parameters["couplings"] *= 3

with tempfile.TemporaryDirectory() as directory:
    checkpoint = os.path.join(directory, "ising.h5")
    with hdf5.archive(checkpoint, "w") as archive:
        simulation.save(archive)
        archive["metadata/matrix"] = np.ones((2, 3), dtype="f8")

    restored = ising_c.sim(parameters)
    with hdf5.archive(checkpoint, "r") as archive:
        restored.load(archive)
        np.testing.assert_array_equal(archive["metadata/matrix"], np.ones((2, 3)))

    after = restored.collectResults()
    np.testing.assert_array_equal(restored.parameters["couplings"], [3., 6.])
    assert restored.parameters["label"] == "checkpoint"
    changed = restored.parameters["couplings"]
    changed[0] = 4
    assert restored.parameters["couplings"][0] == 3
    restored.parameters["couplings"] = changed
    assert restored.parameters["couplings"][0] == 4
    assert restored.resultNames() == simulation.resultNames()
    assert after["Magnetization"].count == before["Magnetization"].count
    np.testing.assert_array_equal(after["Magnetization"].mean, before["Magnetization"].mean)
    np.testing.assert_array_equal(after["Magnetization"].batch_sums, before["Magnetization"].batch_sums)
    np.testing.assert_array_equal(after["Magnetization"].batch_counts, before["Magnetization"].batch_counts)

    target = ising_c.sim(ngs.params({"SEED": 19, "SWEEPS": 3, "label": "kept"}))
    assert target.run(lambda: False)
    handle = target.measurements["Magnetization"]
    saved_result = target.collectResults()["Magnetization"]
    random = copy.deepcopy(target.random)
    random_values = [random() for _ in range(17)]
    broken = os.path.join(directory, "invalid-app.h5")
    for fault in ("state", "sweeps", "count", "missing"):
        shutil.copyfile(checkpoint, broken)
        with h5py.File(broken, "a") as archive:
            if fault == "state":
                archive["checkpoint/state"][()] = 0.
            elif fault == "sweeps":
                archive["checkpoint/sweeps"][()] = 11
            elif fault == "count":
                archive["measurements/Magnetization/batch/count"][0] += 1
            else:
                del archive["checkpoint/state"]
        with hdf5.archive(broken, "r") as archive:
            try:
                target.load(archive)
            except (ValueError, hdf5.ArchiveError):
                pass
            else:
                raise AssertionError(f"invalid {fault} checkpoint accepted")
            assert archive.context == "/"
        assert target.parameters["label"] == "kept" and target.parameters["SWEEPS"] == 3
        assert target.measurements["Magnetization"] is handle
        current = target.collectResults()["Magnetization"]
        assert current.count == saved_result.count == 3
        np.testing.assert_array_equal(current.batch_sums, saved_result.batch_sums)
        np.testing.assert_array_equal(current.batch_counts, saved_result.batch_counts)
        random = copy.deepcopy(target.random)
        assert [random() for _ in range(17)] == random_values

    # `archive[path] = simulation` must write exactly what simulation.save()
    # writes. It reaches save() through the archive-savable marker that this
    # class inherits from mcbase -- but mcbase's own bound save() is a
    # deliberately *non-virtual* qualified call (so that a Python subclass
    # calling super().save(ar) does not re-enter its own override). If the
    # exporter ever stopped binding save() on the derived type, that inherited
    # non-virtual base save is what would run, and this spelling would silently
    # checkpoint the base state only -- a partial write where the previous
    # behaviour was a loud "Unsupported type". Compare the two trees.
    def entries(archive, path="/", found=None, depth=0):
        found = [] if found is None else found
        if depth > 12:
            return found
        for child in archive.list_children(path):
            child_path = path.rstrip("/") + "/" + child
            found.append(child_path)
            if archive.is_group(child_path):
                entries(archive, child_path, found, depth + 1)
        return found

    trees = {}
    for label, write in (("explicit", lambda sim, ar: sim.save(ar)),
                         ("setitem", lambda sim, ar: ar.__setitem__("/", sim))):
        fresh = ising_c.sim(parameters)
        assert fresh.run(lambda: False)
        target = os.path.join(directory, label + ".h5")
        with hdf5.archive(target, "w") as archive:
            write(fresh, archive)
        with hdf5.archive(target, "r") as archive:
            trees[label] = sorted(entries(archive))

    assert trees["setitem"] == trees["explicit"], (
        "archive['/'] = simulation wrote a different tree than "
        "simulation.save(archive):\n"
        f"  only in explicit: {sorted(set(trees['explicit']) - set(trees['setitem']))}\n"
        f"  only in setitem:  {sorted(set(trees['setitem']) - set(trees['explicit']))}")
    # Guard against both spellings degenerating to base-only state.
    assert any(entry.endswith("/checkpoint/sweeps") for entry in trees["setitem"]), \
        trees["setitem"]

faulthandler.cancel_dump_traceback_later()
print("downstream nanobind export: ok")
