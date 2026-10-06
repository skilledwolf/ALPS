# Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
"""Serial/MPI native chains must retain identical evidence and continuation."""
import os
import subprocess
import shlex

import h5py
import numpy as np
import pytest
from pyalps.run_io import write_run_file
from conftest import alps_program


@pytest.fixture
def launcher():
    value = os.environ.get("ALPS_MPIEXEC")
    if not value:
        pytest.skip("set ALPS_MPIEXEC to test an MPI-enabled application build")
    return value


def invoke(launcher, executable, run, processes=1, success=True):
    command = [executable, str(run)]
    if processes > 1:
        command = [launcher, *shlex.split(os.environ.get("ALPS_MPIEXEC_ARGS", "")), "-n", str(processes), *command]
    result = subprocess.run(command, capture_output=True, text=True, timeout=60)
    assert (result.returncode == 0) == success, result.stdout + result.stderr
    return result


def compare(left, right):
    with h5py.File(left) as a, h5py.File(right) as b:
        def check(name, obj):
            assert name in b["simulation"]
            other = b["simulation"][name]
            assert type(obj) is type(other)
            for key in obj.attrs:
                np.testing.assert_array_equal(obj.attrs[key], other.attrs[key], err_msg=name + "/@" + key)
            if isinstance(obj, h5py.Dataset):
                np.testing.assert_array_equal(obj[()], other[()], err_msg=name)
        a["simulation"].visititems(check)


@pytest.mark.parametrize("app", ["simplemc", "spinmc", "qwl", "worm", "dirloop_sse", "loop"])
@pytest.mark.parametrize("rng", ["mt19937", "lagged_fibonacci607"])
@pytest.mark.parametrize("chains", [1, 3])
def test_mpi_native_chains_and_cross_process_restart(launcher, tmp_path, app, rng, chains):
    executable = alps_program(app)
    p = dict(LATTICE="chain lattice", L=5, T=1.8, SWEEPS=37, THERMALIZATION=3)
    if app == "qwl":
        p = dict(LATTICE="chain lattice", L=4, J=1., CUTOFF=12, SWEEPS=3000,
                 NUMBER_OF_WANG_LANDAU_STEPS=3)
    elif app in ("worm", "dirloop_sse", "loop"):
        p = dict(LATTICE="chain lattice", MODEL="spin", L=4, J=1., T=1.,
                 SWEEPS=300, THERMALIZATION=100, SKIP=3)
        if app == "dirloop_sse":
            p.update(INITIAL_CUTOFF=64, **{"MEASURE[Green Function]":True})
    else:
        p.update(ALGORITHM="xy") if app == "simplemc" else p.update(MODEL="O(4)", UPDATE="cluster")
    def run(name, **options):
        execution = dict(seed=2873, bins=8, chains=chains, rng=rng)
        execution.update(options.pop("execution", {}))
        return write_run_file(tmp_path / (name + ".toml"), parameters=p,
            execution=execution, output={"results": name + ".h5", "checkpoint": name + ".checkpoint.h5"}, **options)
    serial, mpi = run("serial"), run("mpi", execution={"checkpoint_interval": 1e-12})
    invoke(launcher, executable, serial)
    invoke(launcher, executable, mpi, processes=2)
    compare(tmp_path / "serial.h5", tmp_path / "mpi.h5")
    compare(tmp_path / "serial.checkpoint.h5", tmp_path / "mpi.checkpoint.h5")
    stopped = run("stopped", execution={"max_sweeps": 13})
    invoke(launcher, executable, stopped, processes=2)
    resumed = run("resumed", input={"checkpoint": "stopped.checkpoint.h5"})
    invoke(launcher, executable, resumed)
    compare(tmp_path / "serial.h5", tmp_path / "resumed.h5")
    repartitioned = run("repartitioned", input={"checkpoint": "stopped.checkpoint.h5"})
    invoke(launcher, executable, repartitioned, processes=3)
    compare(tmp_path / "serial.h5", tmp_path / "repartitioned.h5")


def test_mpi_nonroot_failure_preserves_scientific_output(launcher, tmp_path):
    executable = alps_program("simplemc")
    result = tmp_path / "result.h5"
    checkpoint = tmp_path / "checkpoint.h5"
    result.write_bytes(b"existing result")
    checkpoint.write_bytes(b"existing checkpoint")
    collision = tmp_path / "snapshot.clone2.1.vtk"
    collision.write_bytes(b"existing snapshot")
    run = write_run_file(tmp_path / "failure.toml", parameters=dict(ALGORITHM="ising", LATTICE="chain lattice", L=5, T=2., SWEEPS=37, THERMALIZATION=3),
        execution=dict(seed=31, chains=2, bins=8, snapshot_interval=1),
        output=dict(results="result.h5", checkpoint="checkpoint.h5", snapshot_prefix="snapshot"))
    error = invoke(launcher, executable, run, processes=2, success=False)
    assert collision.name in error.stderr
    assert result.read_bytes() == b"existing result"
    assert checkpoint.read_bytes() == b"existing checkpoint"
    assert collision.read_bytes() == b"existing snapshot"
