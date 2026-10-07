# Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
"""The native loop covers the physical and disorder model of the former loop_single."""
from itertools import product
import subprocess
import xml.etree.ElementTree as ET

import h5py
import numpy as np
import pytest

from conftest import alps_program
from pyalps.run_io import execute, write_run_file, write_run_files
from test_native_mpi import compare


def graph(parameters):
    command = [alps_program('lattice2xml'), parameters['LATTICE']]
    command += [f'{key}={parameters[key]}'
                for key in ('L', 'DEPLETION', 'DEPLETION_SEED') if key in parameters]
    result = subprocess.run(command, text=True, capture_output=True, timeout=30)
    assert result.returncode == 0, result.stderr
    xml = ET.fromstring(result.stdout)
    signs = [(-1.) ** sum(map(int, vertex.findtext('COORDINATE').split()))
             for vertex in xml.findall('VERTEX')]
    bonds = [(int(edge.attrib['source']) - 1, int(edge.attrib['target']) - 1)
             for edge in xml.findall('EDGE')]
    return signs, bonds


def exact_observables(signs, bonds, volume, beta):
    # H = sum_b S_i.S_j. Parallel periodic bonds are separate terms.
    states = np.array(list(product((-.5, .5), repeat=len(signs))))
    indices = {tuple(state): i for i, state in enumerate(states)}
    h = np.zeros((len(states), len(states)))
    for i, state in enumerate(states):
        for left, right in bonds:
            h[i, i] += state[left] * state[right]
            if state[left] != state[right]:
                flipped = state.copy()
                flipped[left], flipped[right] = flipped[right], flipped[left]
                h[i, indices[tuple(flipped)]] += .5
    energies, vectors = np.linalg.eigh(h)
    probabilities = np.exp(-beta * (energies - energies.min()))
    probabilities /= probabilities.sum()
    uniform = states.sum(axis=1)
    staggered = states @ signs
    equal_time = lambda values: (values @ (vectors**2)) @ probabilities

    # The staggered magnetization does not commute with H. Its static
    # susceptibility is the imaginary-time Kubo integral, not beta <M_s^2>.
    matrix = vectors.T @ (staggered[:, None] * vectors)
    gaps = energies[:, None] - energies[None, :]
    degenerate = np.abs(gaps) < 1e-12
    kernel = np.divide(probabilities[None, :] - probabilities[:, None], gaps,
                       out=np.zeros_like(gaps), where=~degenerate)
    kernel[degenerate] = np.broadcast_to(beta * probabilities[:, None], gaps.shape)[degenerate]
    energy = energies @ probabilities
    return {'Volume': volume, 'Energy': energy, 'Energy Density': energy / volume,
            'Staggered Magnetization^2': equal_time(staggered**2),
            'Susceptibility': beta * equal_time(uniform**2) / volume,
            'Staggered Susceptibility': np.sum(matrix**2 * kernel) / volume}


def sections(name, parameters, *, disorder_seed=7, budget=0, checkpoint=None, chains=3):
    return dict(parameters=dict(ALGORITHM='loop', MODEL='spin', local_S=.5,
                                J=1., T=.8, THERMALIZATION=1000, SWEEPS=40000)
                           | parameters,
                execution=dict(seed=137, disorder_seed=disorder_seed,
                               chains=chains, bins=32, max_sweeps=budget),
                input=dict(checkpoint=checkpoint) if checkpoint else None,
                output=dict(results=name + '.h5', checkpoint=name + '.checkpoint.h5'))


def parameters_in(group):
    return {entry['name'].asstr()[()]: entry['value'][()]
            for entry in group['entries'].values()}


@pytest.mark.parametrize('parameters', [
    dict(LATTICE='chain lattice', L=4),
    dict(LATTICE='depleted square lattice', L=2, DEPLETION=.2, DEPLETION_SEED=7),
])
def test_loop_single_observables_match_exact_diagonalization(tmp_path, parameters):
    signs, bonds = graph(parameters)
    expected = exact_observables(signs, bonds, volume=4., beta=1 / .8)
    run = write_run_file(tmp_path / 'physical.toml', **sections('physical', parameters, chains=1))
    execute(alps_program('loop'), run)
    with h5py.File(tmp_path / 'physical.h5') as ar:
        results = ar['simulation/results']
        assert results['Number of Sites/mean/value'][0] == len(signs)
        for name, exact in expected.items():
            mean = results[name + '/mean/value'][0]
            error = results[name + '/mean/error'][0]
            assert abs(mean - exact) < max(.02, 6 * error), (name, mean, error, exact)
    if 'DEPLETION' in parameters:
        # Density and susceptibilities use the undepleted geometric volume,
        # while the Hilbert space contains only the surviving physical sites.
        assert len(signs) == 3
        assert expected['Staggered Susceptibility'] != pytest.approx(
            expected['Staggered Magnetization^2'] / (.8 * 4.))


def test_disorder_realizations_are_explicit_jobs_with_exact_restart(tmp_path):
    # Old NUM_CLONES changed the quenched disorder between clones. Native
    # execution.chains are independent MC histories of one fixed model, so
    # different disorder realizations are explicit jobs with recorded seeds.
    seeds = [7, 17, 27]
    def jobs(prefix, seeds, budget=0, resume=False):
        return write_run_files(tmp_path / prefix, [
            sections(prefix + str(seed),
                     dict(LATTICE='depleted square lattice', L=2, DEPLETION=.2,
                          DEPLETION_SEED=seed, THERMALIZATION=17, SWEEPS=113),
                     disorder_seed=seed, budget=budget, chains=2,
                     checkpoint=f'partial{seed}.checkpoint.h5' if resume else None)
            for seed in seeds])
    executable = alps_program('loop')
    execute(executable, jobs('full', seeds))
    # One realization suffices to show that a restart keeps its disorder.
    execute(executable, jobs('partial', seeds[:1], budget=31))
    execute(executable, jobs('resumed', seeds[:1], resume=True))
    for suffix in ('.h5', '.checkpoint.h5'):
        compare(tmp_path / ('full7' + suffix), tmp_path / ('resumed7' + suffix))
    sites = []
    for seed in seeds:
        with h5py.File(tmp_path / f'full{seed}.h5') as ar:
            assert parameters_in(ar['run_config/execution'])['disorder_seed'] == seed
            clones = ar['simulation/realizations/0/clones']
            assert len(clones) == 2
            volume = ar['simulation/results/Volume/mean/value'][0]
            sites.append(ar['simulation/results/Number of Sites/mean/value'][0])
            assert volume == 4
            for clone in clones.values():
                parameters = parameters_in(clone['sampling_parameters'])
                assert parameters['DISORDER_SEED'] == seed
                assert parameters['DEPLETION_SEED'] == seed
                assert clone['results/Number of Sites/mean/value'][0] == sites[-1]
    assert sites == [3, 4, 3]

