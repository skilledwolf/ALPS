# Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
"""Native sparse and full diagonalization against a dense numpy reference."""
from functools import reduce
import subprocess

import numpy as np
import pytest
import pyalps
from pyalps.run_io import execute, write_run_file
from conftest import alps_program, launch

L = 4
PARAMETERS = dict(LATTICE="chain lattice", MODEL="spin", local_S=.5, J=1, L=L, h=.25, CONSERVED_QUANTUMNUMBERS="Sz")


def site(operator, i):
    return reduce(np.kron, [operator if k == i else np.eye(2) for k in range(L)])


SZ = [site(np.diag([.5, -.5]), i) for i in range(L)]
SPLUS = [site(np.array([[0., 1.], [0., 0.]]), i) for i in range(L)]
HAMILTONIAN = sum(SZ[i] @ SZ[(i + 1) % L] + .5 * (SPLUS[i] @ SPLUS[(i + 1) % L].T + SPLUS[i].T @ SPLUS[(i + 1) % L])
                  for i in range(L)) - .25 * sum(SZ)
TOTAL = np.diag(sum(SZ))


def block(sz):
    """Levels and eigenvectors of the Sz sector, embedded in the full space."""
    indices = np.flatnonzero(TOTAL == sz)
    levels, vectors = np.linalg.eigh(HAMILTONIAN[np.ix_(indices, indices)])
    states = np.zeros((len(TOTAL), len(levels)))
    states[indices] = vectors
    return levels, states


def run(directory, name, **parameters):
    return write_run_file(directory / (name + ".toml"), parameters=dict(PARAMETERS, **parameters),
                          output=dict(results=name + ".out.h5"))


def measurements(result, sz):
    """The first level's measurements in the Sz sector, at k = 0 with momenta."""
    return {s.props["observable"]: np.asarray(s.y)[0]
            for s in pyalps.flatten(pyalps.loadEigenstateMeasurements([result]))
            if s.props["Sz"] == sz and s.props.get("TOTAL_MOMENTUM", 0) == 0}


def test_spectra_measurements_and_thermodynamics(tmp_path):
    # Momentum sectors use complex Bloch states; without translation symmetry
    # the sectors are real and local measurements are available.
    full, = execute(alps_program("fulldiag"), run(tmp_path, "full", **{
        "MEASURE_CORRELATIONS[Diagonal spin correlations]": "Sz"}))
    launch(alps_program("sparsediag"), run(tmp_path, "sparse", TRANSLATION_SYMMETRY=False, **{
        "MEASURE_AVERAGE[Magnetization]": "Sz", "MEASURE_LOCAL[Local magnetization]": "Sz",
        "MEASURE_STRUCTURE_FACTOR[Structure Factor S]": "Sz"}))
    sparse = str(tmp_path / "sparse.out.h5")
    spectrum = pyalps.loadSpectra([full])[0]
    assert all("TOTAL_MOMENTUM" in sector.props for sector in spectrum)
    np.testing.assert_allclose(np.sort(np.concatenate([sector.y for sector in spectrum])),
                               np.linalg.eigvalsh(HAMILTONIAN), atol=1e-12)
    lowest = {sector.props["Sz"]: sector.y for sector in pyalps.loadSpectra([sparse])[0]}
    assert sorted(lowest) == sorted(np.unique(TOTAL))
    for sz, energies in lowest.items():
        np.testing.assert_allclose(energies, block(sz)[0][:1], atol=1e-10)

    # The singlet ground state (at k = 0) and the lowest Sz = 1 state.
    levels, states = block(0)
    ground = states[:, 0]
    correlations = [ground @ SZ[0] @ SZ[d] @ ground for d in range(L)]
    factor = [sum(np.cos(q * (i - j)) * (ground @ SZ[i] @ SZ[j] @ ground) for i in range(L) for j in range(L)) / L
              for q in 2 * np.pi * np.arange(L) / L]
    measured = measurements(full, 0)
    np.testing.assert_allclose(measured["Energy"], levels[0], atol=1e-12)
    np.testing.assert_allclose(measured["Diagonal spin correlations"], correlations, atol=1e-12)
    measured = measurements(sparse, 0)
    np.testing.assert_allclose(measured["Energy"], levels[0], atol=1e-10)
    np.testing.assert_allclose(measured["Magnetization"], 0, atol=1e-10)
    np.testing.assert_allclose(measured["Structure Factor S"], factor, atol=1e-10)
    excited = block(1)[1][:, 0]
    np.testing.assert_allclose(measurements(sparse, 1)["Local magnetization"],
                               [excited @ SZ[i] @ excited for i in range(L)], atol=1e-10)

    # The spectrum computed at h = 1/4 is shifted to h = 1/2 by each level's Sz.
    plots = pyalps.evaluateFulldiagVersusT(full, alps_program("fulldiag_evaluate"),
                                           DELTA_T=.5, T_MIN=.5, T_MAX=1, H=.5)[0]
    curves = {plot.props["ylabel"]: plot for plot in plots}
    temperatures = np.array([.5, 1.])
    sz = np.concatenate([[value] * np.count_nonzero(TOTAL == value) for value in np.unique(TOTAL)])
    shifted = np.concatenate([block(value)[0] for value in np.unique(TOTAL)]) - .25 * sz
    weights = np.exp(-np.outer(1 / temperatures, shifted - shifted.min()))
    z = weights.sum(axis=1)
    energy, moment = weights @ shifted / z, weights @ sz / z
    expected = {"Energy Density": energy / L,
                "Specific Heat per Site": (weights @ shifted**2 / z - energy**2) / temperatures**2 / L,
                "Magnetization per Site": moment / L,
                "Uniform Susceptibility per Site": (weights @ sz**2 / z - moment**2) / temperatures / L}
    for name, values in expected.items():
        np.testing.assert_array_equal(curves[name].x, temperatures)
        np.testing.assert_allclose(curves[name].y, values, atol=1e-10)


@pytest.mark.parametrize("parameters, message", [
    (dict(LATTICE_LIBRARY="lattices.xml"), "Retired parameter LATTICE_LIBRARY"),
    (dict(GRAPH="square lattice"), "exactly one parameters.LATTICE or parameters.GRAPH"),
    (dict(NUMBER_EIGENVALUES=0), "NUMBER_EIGENVALUES"),
])
def test_invalid_runs_are_rejected(tmp_path, parameters, message):
    result = subprocess.run([alps_program("sparsediag"), "--validate", run(tmp_path, "run", **parameters)],
                            capture_output=True, text=True, timeout=60)
    assert result.returncode == 1 and message in result.stderr, result.stderr


def test_runs_cannot_share_results(tmp_path):
    first = write_run_file(tmp_path / "first.toml", parameters=PARAMETERS, output=dict(results="same.h5"))
    second = write_run_file(tmp_path / "second.toml", parameters=PARAMETERS, output=dict(results="same.h5"))
    result = subprocess.run([alps_program("fulldiag"), "--validate", first, second],
                            capture_output=True, text=True, timeout=60)
    assert result.returncode == 1 and "Output paths must be distinct" in result.stderr, result.stderr
