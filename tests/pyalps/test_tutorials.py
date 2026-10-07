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
from conftest import alps_program

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
    executable = alps_program("dmft")
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


QMC_SCRIPTS = scripts(TUTORIALS / '03-mc', "execute('worm'", "execute('dirloop_sse'", "execute('loop'") + scripts(
    TUTORIALS / '11-notebook' / 'ja', "execute('worm'", "execute('dirloop_sse'", "execute('loop'")
QMC_NOTEBOOKS = sorted((TUTORIALS / '11-notebook').rglob('MC-0[23458]*.ipynb'))
# Validating 8d's 32 runs on lattices up to 128x64 takes minutes on slower
# runners; 8a-8c run the same coupled-ladder loop scans at smaller sizes.
LARGE_SCANS = {TUTORIALS / '03-mc/08-quantum-phase-transition/tutorial8d.py'}


@pytest.mark.parametrize('source', [pytest.param(path, marks=pytest.mark.slow) if path in LARGE_SCANS else path
                                    for path in QMC_SCRIPTS + QMC_NOTEBOOKS],
                         ids=lambda path: str(path.relative_to(TUTORIALS)))
def test_native_qmc_tutorial_end_to_end(source, tmp_path, monkeypatch):
    """Validate complete scans, then run representative points through analysis."""
    import tomllib
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import pyalps
    from pyalps import run_io
    original_execute = run_io.execute
    outputs = []

    def execute(app, job):
        executable = alps_program(app)
        _, runs = run_io.read_job_manifest(job)
        # Validation builds every lattice of the scan; large ones take minutes on slow runners.
        checked = subprocess.run([executable, '--validate', *map(str, runs)],
                                 capture_output=True, text=True, timeout=600)
        assert checked.returncode == 0, checked.stdout + checked.stderr
        selected = [runs[i] for i in sorted({0, len(runs)//2, len(runs)-1})]
        # The gap fit needs at least three temperatures in each coupling group.
        documents = {path: tomllib.loads(path.read_text()) for path in runs}
        if app == 'loop' and all('T' in d['parameters'] and 'J2' in d['parameters']
                                 for d in documents.values()):
            groups = {}
            for path, document in documents.items():
                groups.setdefault(document['parameters']['J2'], []).append(path)
            selected = [path for group in groups.values() for path in group[:3]]
        for path in selected:
            document = tomllib.loads(path.read_text())
            p = document['parameters']
            p.update(L=4, THERMALIZATION=1000, SWEEPS=2000)
            if 'W' in p:
                p['W'] = 2
            if 'BETA' in p:
                p['BETA'] = 8.
            document.setdefault('execution', {})['bins'] = 16
            run_io.write_run_file(path, overwrite=True,
                **{key:document[key] for key in ('parameters', 'input', 'output', 'execution') if key in document})
        files = original_execute(executable, selected)
        outputs.extend(files)
        return files

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(run_io, 'execute', execute)
    monkeypatch.setattr(plt, 'show', lambda: None)
    try:
        if source.suffix == '.py':
            runpy.run_path(str(source), run_name='__main__')
        else:
            namespace = {'__name__':'__main__'}
            for code in code_cells(source):
                exec(compile(code, str(source), 'exec'), namespace)
        assert outputs
        for data in pyalps.flatten(pyalps.loadMeasurements(outputs, 'Energy')):
            assert np.isfinite(data.y[0].mean)
            assert np.isfinite(data.y[0].error)
    finally:
        plt.close('all')
