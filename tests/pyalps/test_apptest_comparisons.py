"""Comparison contracts across the application loaders' different layouts."""
import copy

import numpy as np
import pytest
import pyalps
from pyalps import apptest
from pyalps.cxx.pymcdata_c import MCScalarData


def dataset(values, name="energy", filename="test.h5"):
    return pyalps.DataSet(y=values, props={"observable": name, "filename": filename})


def loaders(monkeypatch, test, ref, eigen=None):
    def measurements(files):
        return copy.deepcopy(test if files == "test" else ref)

    def eigenstates(files):
        if eigen is None:
            raise RuntimeError("no eigenstate measurements")
        return copy.deepcopy(eigen[0] if files == "test" else eigen[1])

    monkeypatch.setattr(pyalps, "loadMeasurements", measurements)
    monkeypatch.setattr(pyalps, "loadEigenstateMeasurements", eigenstates)


@pytest.mark.parametrize("mixed", [False, True])
def test_mc_scalar_and_vector_tolerances(monkeypatch, mixed):
    # The worst absolute difference selects the reported vector tolerance.
    test_values = [MCScalarData(1.0, 0.03), MCScalarData(2.0, 0.05)]
    ref_values = [MCScalarData(1.1, 0.04), MCScalarData(2.4, 0.12)]
    test = [dataset([test_values[0]], "scalar"),
            dataset(test_values if mixed else [test_values], "vector")]
    ref = [dataset([ref_values[0]], "scalar"),
           dataset(ref_values if mixed else [ref_values], "vector")]
    loaders(monkeypatch, [test], [ref])
    compare = apptest.compareMixed if mixed else apptest.compareMC
    passed, results = compare("test", "ref")
    scalar, vector = results[0]
    assert scalar["tolerance"] == pytest.approx(0.10)
    # Preserve the existing factor-squared MC vector contract in this refactor.
    assert vector["tolerance"] == pytest.approx(0.52)
    assert vector["difference"] == pytest.approx(0.4)
    assert vector["passed"]
    assert passed == (scalar["passed"] and vector["passed"])


def test_mixed_deterministic_values_and_eigenstate_append(monkeypatch):
    test = [[dataset([1.0], "mc_or_scalar")]]
    ref = [[dataset([1.0], "mc_or_scalar")]]
    eigen = ([[dataset([-10.0, 20.0], "vector")]],
             [[dataset([-10.0 + 1e-10, 20.0], "vector")]])
    loaders(monkeypatch, test, ref, eigen)
    passed, results = apptest.compareMixed("test", "ref", whatlist=["vector"])
    assert not passed
    assert len(results[0]) == 1
    # Vector epsilon tolerances historically use the signed reference value.
    assert results[0][0]["tolerance"] == pytest.approx(2e-11)


def test_epsilon_sectors_and_scalar_relative_tolerance(monkeypatch):
    test = [[[dataset([-10.0])], [dataset([np.array([1.0, 2.0])], "vector")]]]
    ref = [[[dataset([-10.0 + 5e-11])], [dataset([np.array([1.0, 2.0])], "vector")]]]
    loaders(monkeypatch, [], [], (test, ref))
    passed, results = apptest.compareEpsilon("test", "ref")
    assert passed
    assert len(results) == 1 and len(results[0]) == 2
    assert results[0][0]["tolerance"] == pytest.approx(1e-10)
    assert results[0][1]["difference"] == 0


def test_missing_observable_and_mismatched_tasks(monkeypatch, capsys):
    loaders(monkeypatch, [[dataset([1.0])]], [[dataset([1.0])]])
    with pytest.raises(SystemExit, match="1"):
        apptest.compareMixed("test", "ref", whatlist=["missing"])
    assert "missing" in capsys.readouterr().out
    loaders(monkeypatch, [[dataset([1.0])]], [])
    with pytest.raises(Exception, match="number of tasks"):
        apptest.compareMC("test", "ref")


def test_empty_eigenstate_data(monkeypatch, capsys):
    loaders(monkeypatch, [], [], ([], []))
    assert apptest.compareEpsilon("test", "ref") is None
    assert "returned an empty list" in capsys.readouterr().out
