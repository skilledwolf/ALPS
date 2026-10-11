"""Prototype acceptance tests against an installed MPI/OpenMP wheel."""
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import xml.etree.ElementTree as ET

import pytest


def run(args, cwd):
    env = os.environ.copy()
    env['OMP_NUM_THREADS'] = '2'
    env['OPENBLAS_NUM_THREADS'] = '1'
    result = subprocess.run(args, cwd=cwd, env=env, capture_output=True, text=True, timeout=180)
    assert result.returncode == 0, result.stdout + result.stderr
    return result.stdout + result.stderr


@pytest.mark.parametrize('ranks,threads', [(1, 1), (1, 2), (2, 1), (2, 2)])
def test_loop_execution_modes(tmp_path, ranks, threads):
    (tmp_path / 'input').write_text('''LATTICE="chain lattice"
MODEL="spin"
ALGORITHM="loop"
local_S=0.5
L=8
J=1
T=1
THERMALIZATION=8
SWEEPS=32
SEED=42
{T=1;}
{T=2;}
''')
    run(['parameter2xml', 'input'], tmp_path)
    args = ['loop', '--write-xml', '--check-interval', '1', '--total-threads', str(ranks * threads)]
    if ranks > 1:
        args = ['mpiexec', '-n', str(ranks)] + args + ['--mpi', '--threads-per-clone', str(threads)]
    output = run(args + ['input.in.xml'], tmp_path)
    print(output)
    if ranks == 1:
        assert re.search(r'total number of thread\(s\)\s*=\s*' + str(threads), output)
    else:
        assert re.search(r'total number of process\(es\)/thread\(s\)\s*=\s*' + str(ranks) + '/' + str(ranks * threads), output)
    results = list(tmp_path.glob('input.task*.out.xml'))
    assert len(results) == 2
    for path in results:
        root = ET.parse(path).getroot()
        values = root.findall('.//SCALAR_AVERAGE/MEAN')
        assert values, path.read_text()
        assert all(math.isfinite(float(value.text)) for value in values)


def test_python_two_ranks(tmp_path):
    output = run(['mpiexec', '-n', '2', sys.executable, '-c',
                  'from mpi4py import MPI; import pyalps.mpi as m; '
                  'assert m.size == 2; assert MPI.COMM_WORLD.allreduce(m.rank) == 1'], tmp_path)
    print(output)


def test_solver_with_caller_mpi(tmp_path):
    source = Path(__file__).with_name('test_binding_surface.py')
    run([sys.executable, '-m', 'pytest', '-v', str(source) + '::test_ctqmc_solvers_restore_python_signal_handlers'], tmp_path)
