# Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
"""Replica exchange with spatially decomposed walkers: layouts, restart, physics."""
import itertools
from pathlib import Path

import numpy as np
import pyalps
import pytest
from pyalps.run_io import write_run_file
from conftest import tutorials_build, variants
from test_native_mpi import compare, invoke, launcher

BETAS = [.2, .5, .9]
EXAMPLE = Path(__file__).resolve().parents[2]/'tutorials/00-examples/mc/exchange'


@pytest.fixture
def executable():
    return str(tutorials_build()/'00-examples/mc/exchange/exchange')


def run(executable, directory, name, *, processes=1, team=1, launcher=None, budget=0,
        checkpoint=None, parallel=None, length=12, sweeps=97, coupling=.7, rng='mt19937',
        lattice='chain lattice', success=True):
    parameters = dict(ALGORITHM='multiple parallel ising; exchange', LATTICE=lattice, L=length,
                      J=coupling, INVERSE_TEMPERATURE_SET=BETAS, SWEEPS=sweeps,
                      THERMALIZATION=13, RANDOM_EXCHANGE=True)
    execution = dict(seed=4111, bins=8, rng=rng, max_sweeps=budget, processes_per_walker=team,
                     parallel=parallel or ('replicas' if processes > 1 else 'chains'))
    path = write_run_file(directory/(name+'.toml'), parameters=parameters, execution=execution,
                          input={'checkpoint': str(checkpoint)} if checkpoint else {},
                          output=dict(results=name+'.h5', checkpoint=name+'.checkpoint.h5'))
    result = invoke(launcher, executable, path, processes=processes, success=success)
    if not success:
        return result
    return directory/(name+'.h5'), directory/(name+'.checkpoint.h5')


def test_serial_ladder_rejects_teams_and_loop(executable, tmp_path):
    results, _ = run(executable, tmp_path, 'serial')
    assert len(pyalps.loadMeasurements([str(results)])) == len(BETAS)
    failure = run(executable, tmp_path, 'chains', team=2, success=False)
    assert 'processes_per_walker requires' in failure.stdout + failure.stderr
    path = write_run_file(tmp_path/'loop.toml', parameters=dict(ALGORITHM='loop; exchange', L=4),
                          output=dict(results='loop.h5'))
    failure = invoke(None, executable, path, success=False)
    assert 'ALGORITHM' in failure.stdout + failure.stderr


def test_mpi_spatial_example_input(executable, launcher, tmp_path):
    source = EXAMPLE/'spatial.toml'
    (tmp_path/'spatial.toml').write_text(source.read_text())
    result = invoke(launcher, executable, tmp_path/'spatial.toml', processes=2)
    assert result.returncode == 0


# Teams x processes per walker; every layout reproduces the serial ladder.
LAYOUTS = variants((4, 2), (2, 1), (3, 3), (2, 2), (4, 4))


@pytest.mark.parametrize('rng', ['mt19937', 'lagged_fibonacci607'])
@pytest.mark.parametrize('processes,team', LAYOUTS)
def test_mpi_layouts_reproduce_the_serial_ladder(executable, launcher, tmp_path, rng, processes, team):
    full, full_state = run(executable, tmp_path, 'serial', rng=rng)
    result, state = run(executable, tmp_path, 'mpi', rng=rng, processes=processes, team=team,
                        launcher=launcher)
    compare(full, result); compare(full_state, state)


def test_mpi_restart_changes_layout(executable, launcher, tmp_path):
    full, full_state = run(executable, tmp_path, 'full')
    _, part = run(executable, tmp_path, 'part', processes=4, team=2, launcher=launcher, budget=41)
    for name, processes, team in (('teams', 3, 3), ('replicas', 3, 1), ('serial', 1, 1)):
        result, state = run(executable, tmp_path, name, processes=processes, team=team,
                            launcher=launcher, checkpoint=part)
        compare(full, result); compare(full_state, state)
    _, serial_part = run(executable, tmp_path, 'serial-part', budget=63)
    result, state = run(executable, tmp_path, 'from-serial', processes=4, team=2,
                        launcher=launcher, checkpoint=serial_part)
    compare(full, result); compare(full_state, state)


@pytest.mark.parametrize('processes,team,lattice,length', [
    pytest.param(3, 2, 'chain lattice', 12, id='indivisible'),
    pytest.param(4, 4, 'chain lattice', 6, id='sites'),
    pytest.param(2, 2, 'square lattice', 4, id='ring'),
])
def test_mpi_invalid_teams_preserve_outputs(executable, launcher, tmp_path, processes, team,
                                            lattice, length):
    results, state = run(executable, tmp_path, 'kept')
    before = results.read_bytes(), state.read_bytes()
    path = write_run_file(tmp_path/'bad.toml', parameters=dict(
        ALGORITHM='multiple parallel ising; exchange', LATTICE=lattice, L=length,
        INVERSE_TEMPERATURE_SET=BETAS, SWEEPS=97),
        execution=dict(parallel='replicas', processes_per_walker=team),
        output=dict(results='kept.h5', checkpoint='kept.checkpoint.h5'))
    failure = invoke(launcher, executable, path, processes=processes, success=False)
    assert failure.stdout + failure.stderr
    assert before == (results.read_bytes(), state.read_bytes())


@pytest.mark.slow
@pytest.mark.parametrize('coupling', [1., -1.])
def test_mpi_teams_match_exact_ring_thermodynamics(executable, launcher, tmp_path, coupling):
    length = 8
    results, _ = run(executable, tmp_path, 'physics', processes=4, team=2, launcher=launcher,
                     length=length, coupling=coupling, sweeps=40000)
    spins = np.array(list(itertools.product((-1., 1.), repeat=length)))
    energy = -coupling*np.sum(spins*np.roll(spins, 1, axis=1), axis=1)
    magnetization = spins.sum(axis=1)
    observables = {'Energy': energy, 'Magnetization^2': magnetization**2,
                   'Magnetization^4': magnetization**4}
    for replica in pyalps.loadMeasurements([str(results)]):
        beta = next(d.props['BETA'] for d in replica)
        values = {d.props['observable']: d.native_result for d in replica}
        weights = np.exp(-beta*(energy-energy.min()))
        for name, value in observables.items():
            exact = np.average(value, weights=weights)
            np.testing.assert_allclose(values[name].mean, [exact],
                                       atol=max(.02*max(1., abs(exact)), 6*values[name].error[0]))
