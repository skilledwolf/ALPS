# Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
"""Tutorial and notebook runs are accepted by their applications.

Each script runs until its first simulation, which is validated instead of
simulated, so the tutorials stay in step with the applications' run files.
"""
import json
import os
from pathlib import Path
import runpy
import shutil
import subprocess
import sys
import types

import pytest

TUTORIALS = Path(__file__).resolve().parents[2] / "tutorials"
DMFT = TUTORIALS / "05-dmft"
HYBRIDIZATION = TUTORIALS / "06-hybridization"
MARKERS = ("pyalps.run_io", "runDMFT")


class Prepared(Exception):
    """Raised in place of the first simulation once its run was accepted."""


def scripts(directory, *markers):
    return [path for path in sorted(directory.rglob("tutorial*.py"))
            if any(marker in path.read_text(encoding="utf-8") for marker in markers)]


def code_cells(notebook):
    document = json.loads(notebook.read_text(encoding="utf-8"))
    cells = document.get("cells") or [cell for sheet in document.get("worksheets", [])
                                       for cell in sheet["cells"]]
    return ["".join(cell.get("source", cell.get("input", []))) for cell in cells
            if cell["cell_type"] == "code"]


NOTEBOOKS = [path for path in sorted((TUTORIALS / "11-notebook").rglob("DMFT-*.ipynb"))
             if any(marker in cell for cell in code_cells(path) for marker in MARKERS)]


@pytest.fixture
def validate_dmft(monkeypatch):
    executable = os.environ.get("ALPS_DMFT_EXECUTABLE")
    if not executable and os.environ.get("ALPS_DIR"):
        candidate = (Path(os.environ["ALPS_DIR"]).resolve().parents[1] / "bin" /
                     ("dmft.exe" if os.name == "nt" else "dmft"))
        executable = str(candidate) if candidate.is_file() else None
    if not executable:
        pytest.skip("requires ALPS_DMFT_EXECUTABLE or an SDK in ALPS_DIR")
    from pyalps import run_io

    def execute(application, runs, **_):
        assert application == "dmft"
        result = subprocess.run([executable, "--validate", os.fspath(runs)],
                                capture_output=True, text=True, timeout=60)
        assert result.returncode == 0, result.stdout + result.stderr
        raise Prepared

    monkeypatch.setattr(run_io, "execute", execute)


def in_copy(directory, tmp_path, monkeypatch):
    work = tmp_path / directory.name
    shutil.copytree(directory, work)
    monkeypatch.chdir(work)
    return work


@pytest.mark.parametrize("script", scripts(DMFT, *MARKERS), ids=lambda path: path.name)
def test_dmft_tutorial_runs_validate(script, tmp_path, monkeypatch, validate_dmft):
    work = in_copy(script.parent, tmp_path, monkeypatch)
    with pytest.raises(Prepared):
        runpy.run_path(str(work / script.name), run_name="__main__")


@pytest.mark.parametrize("notebook", NOTEBOOKS, ids=lambda path: path.name)
def test_dmft_notebook_runs_validate(notebook, tmp_path, monkeypatch, validate_dmft):
    # A notebook shares its input files with the tutorial of the same number.
    in_copy(next(DMFT.glob(notebook.name[5:7] + "-*")), tmp_path, monkeypatch)
    for cell in code_cells(notebook):
        if any(marker in cell for marker in MARKERS):
            with pytest.raises(Prepared):
                exec(compile(cell, notebook.name, "exec"), {"__name__": "__main__"})


@pytest.mark.parametrize("script", scripts(HYBRIDIZATION, "cthyb"), ids=lambda path: path.name)
def test_hybridization_tutorial_runs_are_prepared(script, tmp_path, monkeypatch):
    cthyb = pytest.importorskip("pyalps.cthyb")
    import pyalps

    world = types.SimpleNamespace(barrier=lambda: None)
    mpi = types.SimpleNamespace(rank=0, size=1, world=world)
    monkeypatch.setitem(sys.modules, "pyalps.mpi", mpi)
    monkeypatch.setattr(pyalps, "mpi", mpi, raising=False)

    def solve(run):
        assert run.application == "cthyb"
        raise Prepared

    monkeypatch.setattr(cthyb, "solve", solve)
    work = in_copy(script.parent, tmp_path, monkeypatch)
    with pytest.raises(Prepared):
        runpy.run_path(str(work / script.name), run_name="__main__")
