# Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
"""DMRG on an open chain against a dense numpy reference."""
from functools import reduce
import subprocess

import numpy as np
import pytest
import pyalps
from pyalps.run_io import execute, write_run_file
from conftest import alps_program

L = 8
J = 1.2345678901
PARAMETERS = {"LATTICE": "open chain lattice", "MODEL": "spin", "local_S": .5, "J": J, "L": L,
              "CONSERVED_QUANTUMNUMBERS": "N,Sz", "Sz_total": 0, "SWEEPS": 2,
              "MEASURE_LOCAL[Local magnetization]": "Sz",
              "MEASURE_CORRELATIONS[Diagonal spin correlations]": "Sz"}


def site(operator, i):
    return reduce(np.kron, [operator if k == i else np.eye(2) for k in range(L)])


SZ = [site(np.diag([.5, -.5]), i) for i in range(L)]
SPLUS = [site(np.array([[0., 1.], [0., 0.]]), i) for i in range(L)]
HAMILTONIAN = J * sum(SZ[i] @ SZ[i + 1] + .5 * (SPLUS[i] @ SPLUS[i + 1].T + SPLUS[i].T @ SPLUS[i + 1])
                      for i in range(L - 1))


def run(directory, name, execution=None, **parameters):
    return write_run_file(directory / (name + ".toml"), parameters=dict(PARAMETERS, **parameters),
                          output={"results": name + ".out.h5"}, execution=execution or {})


def test_energies_measurements_and_iterations(tmp_path):
    # Sixteen states per block make the eight-site chain exact.
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    result, = execute(alps_program("dmrg"), run(tmp_path, "chain", NUMBER_EIGENVALUES=2, STATES=[8, 16, 16, 16],
                                                 execution={"temporary_directory": str(scratch)}))
    assert not any(scratch.iterdir())
    sector = np.flatnonzero(np.diag(sum(SZ)) == 0)
    levels, vectors = np.linalg.eigh(HAMILTONIAN[np.ix_(sector, sector)])
    ground = np.zeros(len(HAMILTONIAN))
    ground[sector] = vectors[:, 0]
    measured = {s.props["observable"]: s for s in pyalps.loadEigenstateMeasurements([result])[0]}
    np.testing.assert_allclose(measured["Energy"].y, levels[:2], atol=1e-10)
    np.testing.assert_allclose(measured["Local magnetization"].y[0], 0, atol=1e-10)
    np.testing.assert_allclose(measured["Diagonal spin correlations"].y[0],
                               [ground @ SZ[i] @ SZ[j] @ ground for i in range(L) for j in range(L)], atol=1e-10)
    history = {s.props["observable"]: s.y for s in pyalps.loadMeasurements([result])[0]}
    assert history["Iteration Energy"][-1] == pytest.approx(levels[0], abs=1e-10)
    assert len(history["Iteration Truncation Error"]) == len(history["Iteration Energy"])


@pytest.mark.parametrize("parameters, execution, message", [
    ({"STATES": [8, 16]}, {}, "2*SWEEPS"),
    ({"MAXSTATES": 16}, {"temporary_directory": "missing"}, "execution.temporary_directory"),
])
def test_invalid_runs_are_rejected(tmp_path, parameters, execution, message):
    result = subprocess.run([alps_program("dmrg"), "--validate", run(tmp_path, "run", execution, **parameters)],
                            capture_output=True, text=True, timeout=60)
    assert result.returncode == 1 and message in result.stderr, result.stderr
