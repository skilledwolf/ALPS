"""Regressions for upstream reports #189--#191 and #197."""
from pathlib import Path

import numpy as np
import pytest

import pyalps
from pyalps import apptest, hdf5, math, ngs


@pytest.mark.parametrize("value", [2, 2., np.array([-2., 0., 2.])])
@pytest.mark.parametrize("name,power", [("sq", 2), ("cb", 3)])
def test_math_powers(value, name, power):
    np.testing.assert_array_equal(getattr(math, name)(value), value ** power)


@pytest.mark.parametrize("name,numpy_name,value", [
    ("asin", "arcsin", .5), ("acos", "arccos", .5), ("atan", "arctan", .5),
    ("asinh", "arcsinh", .5), ("acosh", "arccosh", 2.), ("atanh", "arctanh", .5),
])
def test_math_uses_portable_numpy_names(monkeypatch, name, numpy_name, value):
    expected = getattr(np, numpy_name)(np.array([value]))
    # NumPy 1.26 has only the arc* spellings; exercise that surface on NumPy 2 too.
    monkeypatch.delattr(np, name, raising=False)
    np.testing.assert_allclose(getattr(math, name)(np.array([value])), expected)
    assert getattr(math, name)(value) == pytest.approx(expected[0])


@pytest.mark.parametrize("name", ["sq", "cb", "asin", "acosh"])
def test_math_preserves_object_dispatch(name):
    class Custom:
        pass
    setattr(Custom, name, lambda self: "custom")
    assert getattr(math, name)(Custom()) == "custom"


def test_observable_load_preserves_archive_context_on_success_and_failure(tmp_path):
    observables = ngs.observables()
    observables << ngs.RealObservable("Energy")
    observables["Energy"].append(2.)
    filename = str(tmp_path / "observable.h5")
    with hdf5.archive(filename, "w") as archive:
        archive["/observables"] = observables
        archive["/good/value"] = 7
    with hdf5.archive(filename, "r") as archive:
        archive.set_context("/good")
        observables["Energy"].load(archive, "/observables/Energy")
        assert archive.context == "/good"
        with pytest.raises(hdf5.PathNotFound):
            observables["Energy"].load(archive, "/missing")
        assert archive.context == "/good"
        assert archive["value"] == 7


def input_file(tmp_path, script, references):
    props = {"TESTNAME": "regression", "TOLERANCE": "auto",
             "SAVE_OUT_IF_FAIL": "no", "WRITE_RESULTS_TO_FILE": "no"}
    return apptest.writeTestInputFile(str(script), None, props, references, None,
                                     dirname=str(tmp_path))


@pytest.mark.parametrize("status", [1, 7])
def test_failed_script_cannot_pass_after_writing_outputs(tmp_path, monkeypatch, status):
    monkeypatch.chdir(tmp_path)
    script = tmp_path / "simulation.py"
    script.write_text("from pathlib import Path\nPath('result.h5').touch()\n"
                      f"raise SystemExit({status})\n")
    config = input_file(tmp_path, script, ["result.h5"])
    monkeypatch.setattr(apptest, "compareTest", lambda *a, **kw: pytest.fail("compared failed run"))
    assert apptest.runTest(config) is False
    assert Path.cwd() == tmp_path
    with pytest.raises(RuntimeError, match=f"status {status}"):
        apptest.createTest(str(script), outputs=["result.h5"])
    assert Path.cwd() == tmp_path
    assert not (tmp_path / "apptest_simulation.py").exists()


@pytest.mark.parametrize("create", [False, True])
def test_execution_exception_restores_working_directory(tmp_path, monkeypatch, create):
    monkeypatch.chdir(tmp_path)
    script = tmp_path / "simulation.py"
    script.write_text("")
    config = input_file(tmp_path, script, ["result.h5"])
    def fail(*args):
        raise OSError("cannot execute")
    monkeypatch.setattr(pyalps, "executeCommand", fail)
    with pytest.raises(OSError, match="cannot execute"):
        if create:
            apptest.createTest(str(script), outputs=["result.h5"])
        else:
            apptest.runTest(config)
    assert Path.cwd() == tmp_path


@pytest.mark.parametrize("method", ["auto", "compareEpsilon", "compareMC", "compareMixed"])
@pytest.mark.parametrize("outputs", [[], ["one"], ["one", "two", "extra"]])
def test_comparison_rejects_file_count_mismatch(tmp_path, method, outputs):
    config = input_file(tmp_path, tmp_path / "simulation.py", ["one", "two"])
    with pytest.raises(ValueError, match="number of files"):
        apptest.compareTest(config, outputs, str(tmp_path), "now", method)


@pytest.mark.parametrize("counts", [(1, 2), (2, 1), (0, 1), (1, 0)])
def test_epsilon_rejects_task_count_mismatch(monkeypatch, counts):
    data = iter([[object()] * n for n in counts])
    monkeypatch.setattr(pyalps, "loadEigenstateMeasurements", lambda files: next(data))
    with pytest.raises(ValueError, match="number of tasks"):
        apptest.compareEpsilon(["test"], ["reference"])


def test_successful_script_still_reaches_comparison(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    script = tmp_path / "simulation.py"
    script.write_text("from pathlib import Path\nPath('result.h5').touch()\n")
    config = input_file(tmp_path, script, ["result.h5"])
    compared = []
    def compare(config, outputs, directory, start, **kwargs):
        assert len(outputs) == 1 and Path(outputs[0]).is_file()
        compared.append(outputs[0])
        return True
    monkeypatch.setattr(apptest, "compareTest", compare)
    assert apptest.runTest(config) is True
    assert len(compared) == 1
    assert Path.cwd() == tmp_path
    assert not Path(compared[0]).exists()  # Successful scratch output is removed.
