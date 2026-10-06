# Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
"""convert2xml converts parameter files and rejects Monte Carlo runs."""
import struct
import subprocess

import pytest
from conftest import alps_program

MIGRATION = 'alps-hdf5-convert --alea-results'


def convert(directory, name):
    return subprocess.run([alps_program('convert2xml'), name], cwd=directory,
                          capture_output=True, text=True)


def test_parameter_file_becomes_a_job(tmp_path):
    (tmp_path/'parm').write_text('LATTICE="chain lattice"\nL=4\n{ T=0.5; }\n{ T=1; }\n')
    result = convert(tmp_path, 'parm')
    assert result.returncode == 0, result.stdout + result.stderr
    assert (tmp_path/'parm.in.xml').read_text().count('<TASK') == 2
    task = (tmp_path/'parm.task2.in.xml').read_text()
    assert '<PARAMETER name="T">1</PARAMETER>' in task
    assert '<PARAMETER name="L">4</PARAMETER>' in task


# XDR magic numbers of the former scheduler, simulation and run dumps.
@pytest.mark.parametrize('magic', [1, 2, 3], ids=['scheduler', 'simulation', 'run'])
def test_monte_carlo_dumps_are_rejected(tmp_path, magic):
    (tmp_path/'dump').write_bytes(struct.pack('>i', magic) + bytes(12))
    result = convert(tmp_path, 'dump')
    assert result.returncode != 0
    assert MIGRATION in result.stderr
    assert sorted(path.name for path in tmp_path.iterdir()) == ['dump']


def test_monte_carlo_xml_runs_are_rejected(tmp_path):
    (tmp_path/'run.out.xml').write_text('<?xml version="1.0"?>\n<SIMULATION/>\n')
    result = convert(tmp_path, 'run.out.xml')
    assert result.returncode != 0
    assert MIGRATION in result.stderr
    assert sorted(path.name for path in tmp_path.iterdir()) == ['run.out.xml']
