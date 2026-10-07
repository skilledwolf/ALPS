"""Native result dispatch and the source-tree mean/variance command contracts."""
import contextlib
import io
from pathlib import Path
import runpy
import subprocess
import sys

import h5py
import numpy as np
import pytest

from pyalps import alea, hdf5


TOOLS = Path(__file__).resolve().parents[2] / "src/tools/alea"
FAMILIES = [prefix + family for prefix in ("", "Complex")
            for family in ("Mean", "Variance", "Covariance", "Autocorrelation", "Batch")]
FAMILIES += ["EllipticVariance", "EllipticCovariance"]


def write_result(filename, family="Batch", path="/simulation/results/value"):
    acc = getattr(alea, family + "Accumulator")(2)
    for i in range(33):
        sample = np.array([i, i*i], dtype=float)
        if family.startswith(("Complex", "Elliptic")):
            sample = sample + 1j*sample[::-1]
        acc << sample
    result = acc.result()
    with hdf5.archive(str(filename), "a") as archive:
        result.save(archive, path)
    return result


def command(name, *args, process=False):
    """Run an analysis command like its console script.

    Starting Python and importing pyalps dominates a command, so most checks
    run the script in this interpreter; an escaping exception fails the test
    as a traceback would. process=True starts it as its own process.
    """
    script = TOOLS/(name + ".py")
    if process:
        return subprocess.run([sys.executable, str(script), *map(str, args)],
                              capture_output=True, text=True, timeout=30)
    stdout, stderr, code = io.StringIO(), io.StringIO(), 0
    argv, path = sys.argv, list(sys.path)
    sys.argv, sys.path[:0] = [str(script), *map(str, args)], [str(TOOLS)]
    try:
        with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
            runpy.run_path(str(script), run_name="__main__")
    except SystemExit as exit:
        code = exit.code if isinstance(exit.code, int) else int(exit.code is not None)
    finally:
        sys.argv, sys.path[:] = argv, path
    return subprocess.CompletedProcess(sys.argv, code, stdout.getvalue(), stderr.getvalue())


@pytest.mark.parametrize("family", FAMILIES)
@pytest.mark.parametrize("native", [False, True])
def test_result_dispatch_preserves_evidence(tmp_path, family, native):
    filename = tmp_path / "result.h5"
    expected = write_result(filename, family)
    with hdf5.archive(str(filename)) as ar:
        if native:
            with ar.native() as borrowed:
                result = alea.read_result(borrowed, "/simulation/results/value")
        else:
            result = alea.read_result(ar, "/simulation/results/value")
    assert type(result) is type(expected)
    assert result.count == expected.count
    for field in ("mean", "variance", "covariance", "batch_counts", "batch_sums", "error"):
        if hasattr(expected, field):
            np.testing.assert_array_equal(getattr(result, field), getattr(expected, field))


@pytest.mark.parametrize("name", ["mean", "variance"])
def test_read_only_and_write_options(tmp_path, name):
    filename = tmp_path / "input.h5"
    # Slash in a logical observable name must not be encoded twice when listing.
    from pyalps import hdf5_name_encode
    path = "/custom/" + hdf5_name_encode("Energy/site")
    expected = write_result(filename, "ComplexBatch", path)
    before = filename.read_bytes()
    result = command(name, "-v", "-p", "/custom", filename)
    assert result.returncode == 0, result.stderr
    assert "Energy/site" in result.stdout
    assert filename.read_bytes() == before
    result = command(name, "-w", "-n", "Energy/site", "-p", "/custom", filename)
    assert result.returncode == 0, result.stderr
    with h5py.File(filename) as archive:
        np.testing.assert_array_equal(archive[path + "/" + name + "/value"][()], getattr(expected, name))
    with hdf5.archive(str(filename)) as archive:
        reread = alea.read_result(archive, path)
    np.testing.assert_array_equal(reread.batch_sums, expected.batch_sums)
    np.testing.assert_array_equal(reread.batch_counts, expected.batch_counts)


def test_mean_only_variance_is_not_invented_or_partially_written(tmp_path):
    filename = tmp_path / "mixed.h5"
    write_result(filename, "Batch", "/simulation/results/first")
    write_result(filename, "Mean", "/simulation/results/second")
    result = command("variance", "-w", filename, process=True)
    assert result.returncode != 0
    assert "mean only" in result.stderr
    with h5py.File(filename) as archive:
        assert "variance" not in archive["/simulation/results/first"]
        assert "variance" not in archive["/simulation/results/second"]


@pytest.mark.parametrize("kind", [None, 123, True, [5, 5]])
def test_reject_legacy_checkpoint_or_invalid_kind(tmp_path, kind):
    filename = tmp_path / "invalid.h5"
    write_result(filename)
    with h5py.File(filename, "a") as archive:
        group = archive["simulation/results/value"]
        del group.attrs["kind"]
        if kind is not None:
            group.attrs["kind"] = kind
    result = command("mean", filename, process=kind is None)
    assert result.returncode != 0
    assert "Traceback" not in result.stderr
    if kind is None:
        assert "alps-hdf5-convert" in result.stderr


def test_missing_input_is_not_created_by_write(tmp_path):
    filename = tmp_path / "missing.h5"
    assert command("mean", "-w", filename).returncode != 0
    assert not filename.exists()


def test_commands_read_converted_released_core_results(tmp_path):
    import json
    root = Path(__file__).resolve().parents[2]
    fixture = root / "tests/cli/fixtures/alpscore-v2.3.3-alea.h5"
    selections = json.loads(fixture.with_suffix('.json').read_text())["selections"]
    output = tmp_path / "converted.h5"
    args = [sys.executable, str(root/"src/tools/hdf5/convert.py"), str(fixture), str(output)]
    for name, kind in selections.items():
        args.extend(["--core-alea", kind, "/results/" + name])
    subprocess.run(args, check=True, capture_output=True, text=True, timeout=30)
    for estimate in ("mean", "variance"):
        names = [name for name in selections if not (name == "mean" and estimate == "variance")]
        result = command(estimate, "-v", "-p", "/results", *(f"--name={name}" for name in names), output)
        assert result.returncode == 0, result.stderr
        for name in names:
            assert f" of variable {name} in file " in result.stdout
