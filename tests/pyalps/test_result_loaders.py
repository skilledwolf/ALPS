# Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
"""Contracts for the different result layouts consumed by pyalps.load."""

import numpy as np
import pytest

import pyalps
from pyalps import hdf5
from pyalps.load import Hdf5Loader


@pytest.fixture
def results(tmp_path):
    filename = str(tmp_path / "results.h5")
    with hdf5.archive(filename, "w") as ar:
        ar["/parameters/L"] = 2
        ar["/parameters/T"] = 1.5
        ar["/parameters/MODEL"] = "spin"
        ar["/parameters/couplings"] = [1, 2]
        ar["/simulation/results/Energy/mean/value"] = -2.0
        ar["/simulation/results/Correlations/mean/value"] = [1.0, 0.5]
        ar["/simulation/results/Correlations/labels"] = ["0--0", "0--1"]
        encoded = pyalps.hdf5_name_encode("Density / Site")
        ar[f"/simulation/results/{encoded}/mean/value"] = [0.25, 0.75]
        ar["/simulation/results/Histogram/histogram"] = [3.0, 5.0, 2.0]
        ar["/simulation/results/Histogram/@min"] = -1.0
        ar["/simulation/results/Histogram/@stepsize"] = 0.5
        # Simple-binning arrays retain the existing lattice reshaping convention.
        ar["/simulation/results/Local/mean/value"] = [1.0, 2.0, 3.0, 4.0]
        ar["/simulation/results/Local/mean/error"] = [0.1] * 4
        ar["/spectrum/energies"] = [-2.0, -1.0]
        ar["/quantumnumbers/Sz"] = 0
        for path, shift in (("/spectrum/results", 0),
                            ("/spectrum/sectors/0/results", 10),
                            ("/spectrum/sectors/1/results", 20)):
            ar[path + "/Correlations/mean/value"] = np.arange(6).reshape(2, 3) + shift
            ar[path + "/Correlations/labels"] = ["0", "1", "2"]
            ar[path + "/Energy/mean/value"] = [-2.0 + shift, -1.0 + shift]
        for sector in (0, 1):
            path = f"/spectrum/sectors/{sector}"
            ar[path + "/energies"] = [-2.0 + sector, -1.0 + sector]
            ar[path + "/quantumnumbers/Sz"] = sector
        for iteration in (1, 2):
            path = f"/spectrum/iteration/{iteration}"
            ar[path + "/parameters/T"] = iteration + 2.0
            ar[path + "/results/Energy/mean/value"] = [-2.0 / iteration]
            path = f"/simulation/iteration/{iteration}/results/G_tau"
            ar[path + "/0/mean/value"] = [float(iteration), iteration + 0.5]
            ar[path + "/1/mean/value"] = float(iteration)
        # The last seven binning levels are intentionally excluded by the reader.
        path = "/simulation/results/Energy/timeseries"
        ar[path + "/logbinning"] = np.array([4., 8.] + [0.] * 7)
        ar[path + "/logbinning2"] = np.array([12., 12.] + [0.] * 7)
        ar[path + "/logbinning_counts"] = np.array([4., 4.] + [1.] * 7)
        for step in (1, 2):
            path = f"/timesteps/{step}"
            ar[path + "/parameters/time"] = step * 0.5
            ar[path + "/parameters/T"] = 99.0
            ar[path + "/results/Energy/mean/value"] = -float(step)
    return filename


def test_properties_and_observable_names(results):
    props, = pyalps.loadProperties([results])
    assert props["filename"] == results
    assert props["L"] == 2.0
    assert props["MODEL"] == "spin"
    assert props["couplings"] == [1.0, 2.0]
    names, = pyalps.loadObservableList([results])
    assert "Density / Site" in names
    assert names == props["ObservableList"]


def test_measurement_filtering_order_shapes_and_labels(results):
    requested = ["Density / Site", "Energy", "missing", "Correlations"]
    datasets, = pyalps.loadMeasurements([results], requested)
    assert [d.props["observable"] for d in datasets] == [requested[i] for i in (0, 1, 3)]
    np.testing.assert_array_equal(datasets[0].y, [0.25, 0.75])
    np.testing.assert_array_equal(datasets[1].y, [-2.0])
    np.testing.assert_array_equal(datasets[1].x, [0])
    np.testing.assert_array_equal(datasets[2].x, [0, 1])
    assert datasets[0].props["hdf5_path"].endswith("Density &#47; Site")
    assert all(d.props["filename"] == results for d in datasets)


def test_histogram_and_simple_binning_shapes(results):
    (histogram, local), = pyalps.loadMeasurements([results], ["Histogram", "Local"])
    np.testing.assert_array_equal(histogram.x, [-1.0, -0.5, 0.0])
    np.testing.assert_array_equal(histogram.y, [3.0, 5.0, 2.0])
    np.testing.assert_array_equal(local.y, [[1.0, 2.0], [3.0, 4.0]])
    np.testing.assert_array_equal(local.x, [0, 1, 2, 3])
    assert local.props["origin"] == [0.5, 0.5]
    assert "origin" not in histogram.props


def test_spectra_are_flat_within_each_file(results):
    spectra, = pyalps.loadSpectra([results])
    assert len(spectra) == 3
    assert [d.props["Sz"] for d in spectra] == [0, 0, 1]
    assert [d.props["observable"] for d in spectra] == ["spectrum"] * 3
    np.testing.assert_array_equal(spectra[0].y, [-2.0, -1.0])
    assert isinstance(spectra[0].x, range)
    assert list(spectra[0].x) == [0, 1]


@pytest.mark.parametrize("index", [None, 1])
def test_diagonalization_sector_nesting_and_index(results, index):
    data, = Hdf5Loader().ReadDiagDataFromFile(
        [results], measurements=["Correlations"], index=index
    )
    assert len(data) == 3
    assert isinstance(data[1], list) and isinstance(data[2], list)
    expected = np.arange(6).reshape(2, 3)
    if index is not None:
        expected = expected[index]
    for dataset, shift in ((data[0], 0), (data[1][0], 10), (data[2][0], 20)):
        np.testing.assert_array_equal(dataset.y, expected + shift)
        np.testing.assert_array_equal(dataset.x, [0, 1, 2])
    assert data[2][0].props["Sz"] == 1


def test_iteration_nesting_and_local_parameter_precedence(results):
    data, = pyalps.loadIterationMeasurements([results], "Energy")
    iterations = data[0]
    assert len(iterations) == 2
    for number, (dataset,) in enumerate(iterations, 1):
        assert dataset.props["iteration"] == str(number)
        assert dataset.props["T"] == number + 2.0
        assert dataset.props["L"] == 2.0
        np.testing.assert_array_equal(dataset.y, [-2.0 / number])
    assert data[1][0].props["Sz"] == 0
    assert data[2][0].props["Sz"] == 1


def test_dmft_file_iteration_measurement_nesting(results):
    iterations, = pyalps.loadDMFTIterations([results], measurements=["1", "0"])
    assert len(iterations) == 2
    for number, (scalar, vector) in enumerate(iterations, 1):
        assert scalar.props["observable"] == "1"
        assert scalar.props["iteration"] == str(number)
        np.testing.assert_array_equal(scalar.y, [float(number)])
        np.testing.assert_array_equal(vector.y, [number, number + 0.5])
        np.testing.assert_array_equal(vector.x, [0, 1])


def test_binning_analysis_values(results):
    (dataset,), = pyalps.loadBinningAnalysis([results], "Energy")
    np.testing.assert_allclose(dataset.y, np.sqrt([0.5, 0.5]))
    np.testing.assert_array_equal(dataset.x, [0, 1])
    assert dataset.props["observable"] == "binning analysis of Energy"


def test_time_evolution_flattens_files_and_global_parameters_win(results):
    data = pyalps.load.loadTimeEvolution([results], measurements=["Energy"])
    assert len(data) == 2
    for step, (dataset,) in enumerate(data, 1):
        assert dataset.props["time"] == step * 0.5
        assert dataset.props["T"] == 1.5
        np.testing.assert_array_equal(dataset.y, [-float(step)])


def test_filename_normalization_missing_files_and_file_order(results, tmp_path, capsys):
    data = pyalps.loadMeasurements(
        [str(tmp_path / "missing.h5"), results[:-3] + ".xml", results[:-3]], "Energy"
    )
    assert len(data) == 2
    assert "DOES NOT EXIST" in capsys.readouterr().out
    assert all(group[0].props["filename"] == results for group in data)


def test_file_errors_skip_bad_input_and_continue(results, tmp_path, capsys):
    broken = str(tmp_path / "broken.h5")
    with hdf5.archive(broken, "w") as ar:
        ar["/unrelated"] = 1
    data = pyalps.loadMeasurements([broken, results], "Energy")
    assert len(data) == 1
    assert data[0][0].props["filename"] == results
    assert capsys.readouterr().out


def test_indexed_scalar_with_labels_and_invalid_index(results):
    with hdf5.archive(results, "a") as ar:
        for path in ("/spectrum/results", "/spectrum/sectors/0/results", "/spectrum/sectors/1/results"):
            ar[path + "/Energy/labels"] = ["0"]
    data, = Hdf5Loader().ReadDiagDataFromFile([results], measurements=["Energy"], index=1)
    assert data[0].y == -1.0
    np.testing.assert_array_equal(data[0].x, [0])
    data, = Hdf5Loader().ReadDiagDataFromFile([results], measurements=["Correlations"], index=99)
    assert data[0].y.size == 0
    np.testing.assert_array_equal(data[0].x, [0, 1, 2])


@pytest.mark.parametrize("vector", [False, True])
def test_statistical_measurements_retain_errors_and_counts(results, vector):
    from pyalps import alea

    observable_type = alea.RealVectorTimeSeriesObservable if vector else alea.RealTimeSeriesObservable
    observable = observable_type("Samples")
    samples = np.arange(32.0)
    if vector:
        samples = np.column_stack((samples, 2 * samples))
    for sample in samples:
        observable << sample
    observable.save(results)
    (dataset,), = pyalps.loadMeasurements([results], "Samples")
    data = dataset.y if vector else dataset.y[0]
    assert isinstance(data, alea.MCVectorData if vector else alea.MCScalarData)
    assert data.count == 32
    np.testing.assert_allclose(data.mean, samples.mean(axis=0))
    np.testing.assert_allclose(data.error, observable.error)
    np.testing.assert_array_equal(dataset.x, [0, 1] if vector else [0])


@pytest.mark.parametrize("method", [
    "GetProperties", "ReadSpectrumFromFile", "ReadDiagDataFromFile",
    "ReadBinningAnalysis", "ReadMeasurementFromFile", "ReadDMFTIterations",
])
@pytest.mark.parametrize("fail", [False, True])
def test_archives_close_on_success_and_error(results, monkeypatch, method, fail):
    opened = []
    archive_type = hdf5.archive

    def open_archive(*args, **kwargs):
        archive = archive_type(*args, **kwargs)
        opened.append(archive)
        return archive

    monkeypatch.setattr(hdf5, "archive", open_archive)
    loader = Hdf5Loader()
    if fail:
        def fail_to_read(path):
            raise RuntimeError("unreadable parameters")
        monkeypatch.setattr(loader, "ReadParameters", fail_to_read)
    if fail and method == "ReadDiagDataFromFile":
        with pytest.raises(RuntimeError, match="unreadable parameters"):
            getattr(loader, method)([results])
    else:
        loaded = getattr(loader, method)([results])
        if fail:
            assert loaded == []
        else:
            assert len(loaded) == 1
    assert len(opened) == 1
    assert opened[0].closed


def test_time_evolution_closes_every_archive(results, monkeypatch):
    opened = []
    archive_type = hdf5.archive

    def open_archive(*args, **kwargs):
        archive = archive_type(*args, **kwargs)
        opened.append(archive)
        return archive

    monkeypatch.setattr(hdf5, "archive", open_archive)
    assert len(pyalps.load.loadTimeEvolution([results], measurements=["Energy"])) == 2
    assert opened and all(archive.closed for archive in opened)
