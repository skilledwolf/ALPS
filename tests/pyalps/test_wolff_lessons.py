"""Native Wolff lessons preserve cluster physics and correlated Binder evidence."""
import itertools
import subprocess

import numpy as np
import pytest
import pyalps
from pyalps import alea, hdf5
from conftest import tutorials_build


@pytest.fixture(params=['07-alea', '08-lattice'])
def lesson(request):
    return request.param, str(tutorials_build()/'08-alpsize'/request.param/'wolff')


def config(path, lesson, extra='', output='result.h5', default=False):
    graph = ('LATTICE="square lattice"\n' if default else 'LATTICE="chain lattice"\n') if lesson == '08-lattice' else ''
    parameters = '' if default else 'L=4\nT=2.5\nSWEEPS=100000\nTHERMALIZATION=1000\n'
    path.write_text('[parameters]\n'+graph+parameters+extra+'\n[output]\nresults="'+output+'"\n')


def test_physics_and_joint_binder(lesson, tmp_path):
    name, executable = lesson
    run = tmp_path/'run.toml'
    config(run, name)
    subprocess.run([executable, str(run)], check=True, capture_output=True)
    filename = str(tmp_path/'result.h5')
    with hdf5.archive(filename) as archive:
        joint = alea.read_result(archive, '/simulation/joint')
    assert joint.count == 100000
    # Exact finite systems: square 4x4 (lesson 07), periodic four-site chain (08).
    n = 16 if name == '07-alea' else 4
    spins = np.array(list(itertools.product((-1., 1.), repeat=n)))
    m = spins.mean(axis=1)
    if n == 16:
        grid = spins.reshape(-1, 4, 4)
        energy = -np.sum(grid*(np.roll(grid, 1, 1)+np.roll(grid, 1, 2)), axis=(1, 2))
    else:
        energy = -np.sum(spins*np.roll(spins, 1, 1), axis=1)
    expected = np.average(np.array([m, m*m, m**4]), axis=1, weights=np.exp(-energy/2.5))
    np.testing.assert_allclose(joint.mean, expected, atol=0.015)
    data = {d.props['observable']: d.native_result for d in pyalps.loadMeasurements([filename])[0]}
    ratio = data['Binder Ratio of Magnetization']
    keep = joint.batch_counts > 0
    weights, sums = joint.batch_counts[keep].astype(float), joint.batch_sums[keep]
    total = sums.sum(axis=0)
    mean, count = total/weights.sum(), weights.sum()
    leave = (total-sums)/(count-weights[:, None])
    pseudovalues = count*(mean[1]**2/mean[2])-(count-weights)*(leave[:, 1]**2/leave[:, 2])
    np.testing.assert_allclose(ratio.batch_sums[keep, 0], pseudovalues, rtol=1e-10)
    expected_mean = pseudovalues.sum()/count
    variance = np.sum(weights*(pseudovalues/weights-expected_mean)**2)/(count-np.sum(weights**2)/count)
    np.testing.assert_allclose(ratio.mean, [expected_mean], rtol=1e-12)
    np.testing.assert_allclose(ratio.error, [np.sqrt(variance*np.sum(weights**2)/count**2)], rtol=1e-10)
    for i, observable in enumerate(('Magnetization', 'Magnetization^2', 'Magnetization^4')):
        np.testing.assert_allclose(data[observable].mean, joint.mean[i:i+1], rtol=1e-12)
    diagnostics = pyalps.loadBinningAnalysis([filename])[0]
    assert len(diagnostics) == 3 and all(d.native_result.count == 100000 for d in diagnostics)
    subprocess.run([executable, str(run)], check=True, capture_output=True)
    with hdf5.archive(filename) as archive:
        np.testing.assert_array_equal(alea.read_result(archive, '/simulation/joint').batch_sums, joint.batch_sums)


def test_validation(lesson, tmp_path):
    name, executable = lesson
    run = tmp_path/'run.toml'
    config(run, name)
    subprocess.run([executable, '--validate', str(run)], check=True, capture_output=True)
    assert not (tmp_path/'result.h5').exists()
    schema = subprocess.run([executable, '--schema'], check=True, capture_output=True, text=True)
    assert 'wolff-tutorial' in schema.stdout
    run.write_text(run.read_text().replace('T=2.5', 'T=0.0'))
    assert subprocess.run([executable, '--validate', str(run)], capture_output=True).returncode != 0
    config(run, name, output='run.toml')
    before = run.read_bytes()
    assert subprocess.run([executable, str(run)], capture_output=True).returncode != 0
    assert run.read_bytes() == before


def test_original_square_lesson_means(lesson, tmp_path):
    name, executable = lesson
    run = tmp_path/'run.toml'
    config(run, name, default=True)
    subprocess.run([executable, str(run)], check=True, capture_output=True)
    with hdf5.archive(str(tmp_path/'result.h5')) as archive:
        joint = alea.read_result(archive, '/simulation/joint')
    # Rounded means from each retired wolff.op, not its legacy error estimates.
    expected = ([-0.000368834, 0.626016, 0.408456] if name == "07-alea"
                else [-0.00363082, 0.624276, 0.406649])
    np.testing.assert_allclose(joint.mean, expected, atol=5e-7, rtol=0)


def test_unavailable_binder_keeps_measured_evidence(lesson, tmp_path):
    name, executable = lesson
    run = tmp_path/'run.toml'
    graph = 'LATTICE="chain lattice"\n' if name == '08-lattice' else ''
    length = 4 if name == '08-lattice' else 2
    # Two single-spin flips produce m=.5 then m=0; one leave-out ratio is undefined.
    run.write_text(f'[parameters]\n{graph}L={length}\nT=1e308\nSWEEPS=2\nTHERMALIZATION=0\n'
                   '[execution]\nseed=1\n[output]\nresults="result.h5"\n')
    subprocess.run([executable, str(run)], check=True, capture_output=True)
    with hdf5.archive(str(tmp_path/'result.h5')) as archive:
        joint = alea.read_result(archive, '/simulation/joint')
        assert joint.count == 2
        np.testing.assert_array_equal(joint.mean, [.25, .125, .03125])
        assert not archive.is_group('/simulation/results/Binder Ratio of Magnetization')
        assert 'positive <m^4>' in archive['/simulation/unavailable/Binder Ratio of Magnetization']
