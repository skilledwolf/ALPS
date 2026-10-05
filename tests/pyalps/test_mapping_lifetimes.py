"""Native measurement handles and collected snapshots own their state."""

import copy
import gc

import numpy as np
import pytest

from pyalps import alea, hdf5, ngs


class Simulation(ngs.mcbase):
    def update(self): pass
    def measure(self): pass
    def fraction_completed(self): return 1.0


@pytest.mark.parametrize("remove", ["del", "clear", "replace"])
def test_measurement_handles_survive_collection_mutation_and_owner(tmp_path, remove):
    original = Simulation({"SEED": 1})
    original.measurements["x"] = alea.BatchAccumulator()
    original.measurements["x"] << 1.0 << 3.0
    filename = tmp_path / "state.h5"
    hdf5.save_checkpoint(filename, original.save)
    simulation = Simulation({"SEED": 2})
    with hdf5.archive(filename) as archive:
        simulation.load(archive)
    handle = simulation.measurements["x"]
    handle << 8.0
    np.testing.assert_array_equal(ngs.collectResults(simulation)["x"].mean, [4.0])
    if remove == "del":
        del simulation.measurements["x"]
    elif remove == "clear":
        simulation.measurements.clear()
    else:
        simulation.measurements["x"] = alea.BatchAccumulator()
    del original, simulation
    gc.collect()
    handle << 4.0
    np.testing.assert_array_equal(handle.result().mean, [4.0])
    assert handle.count == 4


def test_collected_results_are_independent_owning_dicts():
    simulation = Simulation({"SEED": 1})
    simulation.measurements["x"] = alea.BatchAccumulator()
    simulation.measurements["x"] << 1.0 << 3.0
    first = ngs.collectResults(simulation)
    assert isinstance(first, dict)
    snapshot = first.pop("x")
    simulation.measurements["x"] << 8.0
    np.testing.assert_array_equal(snapshot.mean, [2.0])
    np.testing.assert_array_equal(ngs.collectResults(simulation)["x"].mean, [4.0])
    del simulation, first
    gc.collect()
    np.testing.assert_array_equal(snapshot.mean, [2.0])


def test_native_accumulator_and_result_copies_are_independent():
    original = alea.BatchAccumulator(2, num_batches=8, base_size=3)
    for i in range(11):
        original << [i, 2*i]
    duplicate = copy.deepcopy(original)
    snapshot = copy.deepcopy(original.result())
    original << [100, 200]
    assert original.count == 12 and duplicate.count == snapshot.count == 11
    np.testing.assert_array_equal(duplicate.result().batch_sums, snapshot.batch_sums)
    duplicate.reset()
    assert snapshot.count == 11


def test_params_iterator_survives_mutation():
    parameters = ngs.params({"first": 1, "second": 2})
    iterator = iter(parameters)
    parameters.clear()
    parameters["new"] = 3
    del parameters
    assert list(iterator) == ["first", "second"]
