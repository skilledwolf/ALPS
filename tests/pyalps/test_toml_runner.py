# Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
"""Explicit jobs use the selected application's native schema and safe argv."""
from pathlib import Path
import subprocess
import threading

import pytest

from pyalps import hdf5, run_config, tools
from pyalps.run_io import execute, write_run_file, write_run_files
from conftest import alps_program


SCHEMA = '''
application="cthyb"
schema_version=1
[parameters.count]
type="int64"
required=true
min=1
[input.data]
type="path"
[output.results]
type="path"
required=true
[output.checkpoint]
type="path"
[execution.seed]
type="int64"
default=42
'''


@pytest.fixture
def commands(monkeypatch):
    recorded = []

    def invoke(arguments, **kwargs):
        assert kwargs["check"] is True
        assert not kwargs.get("shell", False)
        recorded.append(list(arguments))
        return subprocess.CompletedProcess(arguments, 0, stdout=SCHEMA if "--schema" in arguments else "")

    monkeypatch.setattr(tools, "check_existence", lambda _: None)
    monkeypatch.setattr(subprocess, "run", invoke)
    return recorded


def write_job(directory, application="cthyb"):
    return write_run_files(directory / 'batch "quoted"; literal', [
        {"parameters": {"count": 1}, "input": {"data": "data with spaces.h5"},
         "output": {"results": "first result.h5"}},
        {"parameters": {"count": 2}, "output": {"results": "second result.h5"}},
    ], SCHEMA.replace('application="cthyb"', f'application="{application}"'))


def executed(commands):
    return [command for command in commands if "--schema" not in command and "--validate" not in command]


def test_manifest_runs_once_and_returns_actual_outputs_from_another_cwd(tmp_path, monkeypatch, commands):
    directory = tmp_path.resolve() / "runs with spaces"
    manifest = write_job(directory)
    other = tmp_path.resolve() / "working directory"
    other.mkdir()
    monkeypatch.chdir(other)
    application = str(tmp_path / "bin with spaces" / "hybridization")
    results = execute(application, manifest)
    runs = [str(directory / 'batch "quoted"; literal.task1.toml'),
            str(directory / 'batch "quoted"; literal.task2.toml')]
    assert results == [str(directory / "first result.h5"), str(directory / "second result.h5")]
    assert executed(commands) == [[application, run] for run in runs]
    # Every run is validated before the first one starts.
    validations = [i for i, command in enumerate(commands) if "--validate" in command]
    assert sorted(command[-1] for command in commands if "--validate" in command) == runs
    assert max(validations) < commands.index(executed(commands)[0])
    assert Path.cwd() == other
    assert not (other / "literal").exists()


def test_direct_runs_and_mpi_use_only_an_external_launcher(tmp_path, commands):
    tmp_path = tmp_path.resolve()
    paths = [write_run_file(tmp_path / f"run {n}.toml", SCHEMA,
                            parameters={"count": n}, output={"results": f"result {n}.h5"})
             for n in (1, 2)]
    assert execute("hybridization", paths, mpi=2, mpirun="mpi launcher") == [
        str(tmp_path / "result 1.h5"), str(tmp_path / "result 2.h5")]
    assert executed(commands) == [["mpi launcher", "-np", "2", "hybridization", str(path)]
                                 for path in paths]
    assert not any("--mpi" in command for command in commands)


def test_manifest_identity_is_checked_against_selected_executable(tmp_path, commands):
    manifest = write_job(tmp_path, application="ctint")
    with pytest.raises(ValueError, match="does not match executable application"):
        execute("hybridization", manifest)
    assert executed(commands) == []
    assert not any("--validate" in command for command in commands)


def test_only_toml_runs_are_executed(tmp_path, commands):
    with pytest.raises(ValueError, match="TOML"):
        execute("hybridization", [tmp_path / "run.toml", tmp_path / "run.in.xml"])
    manifest = tmp_path / "legacy.job.toml"
    manifest.write_text('[[runs]]\nfile="legacy.in.xml"\n', encoding="utf-8")
    with pytest.raises(ValueError, match="TOML run files"):
        execute("hybridization", manifest)
    with pytest.raises(ValueError, match="mpi"):
        execute("hybridization", tmp_path / "run.toml", mpi=0)
    assert commands == []


def test_repeated_run_aliases_are_not_executed_twice(tmp_path, commands):
    path = write_run_file(tmp_path / "run.toml", SCHEMA,
                          parameters={"count": 1}, output={"results": "result.h5"})
    alias = tmp_path / "alias.toml"
    try:
        alias.symlink_to(path)
    except OSError:
        pytest.skip("file symlinks unavailable")
    with pytest.raises(ValueError, match="duplicate TOML"):
        execute("hybridization", [path, alias])
    assert executed(commands) == []


def test_later_validation_failure_prevents_all_execution(tmp_path, commands, monkeypatch):
    manifest = write_job(tmp_path)
    original = subprocess.run

    def invoke(arguments, **kwargs):
        if "--validate" in arguments and "task2" in arguments[-1]:
            raise subprocess.CalledProcessError(7, arguments)
        return original(arguments, **kwargs)

    monkeypatch.setattr(subprocess, "run", invoke)
    with pytest.raises(subprocess.CalledProcessError) as failure:
        execute("hybridization", manifest)
    assert failure.value.returncode == 7
    assert executed(commands) == []


def test_later_process_failure_propagates_and_stops_remaining_runs(tmp_path, commands, monkeypatch):
    tmp_path = tmp_path.resolve()
    paths = [write_run_file(tmp_path / f"run{n}.toml", SCHEMA,
                            parameters={"count": n}, output={"results": f"result{n}.h5"})
             for n in (1, 2, 3)]
    original = subprocess.run

    def invoke(arguments, **kwargs):
        if len(arguments) == 2 and arguments[-1] == str(paths[1]):
            commands.append(list(arguments))
            raise subprocess.CalledProcessError(9, arguments)
        return original(arguments, **kwargs)

    monkeypatch.setattr(subprocess, "run", invoke)
    with pytest.raises(subprocess.CalledProcessError) as failure:
        execute("hybridization", paths)
    assert failure.value.returncode == 9
    assert executed(commands) == [["hybridization", str(paths[0])], ["hybridization", str(paths[1])]]


def concurrent_runs(directory, count):
    return [write_run_file(directory.resolve() / f"run{n}.toml", SCHEMA,
                           parameters={"count": n}, output={"results": f"result{n}.h5"})
            for n in range(1, count + 1)]


def test_concurrent_runs_overlap_and_return_results_in_input_order(tmp_path, commands, monkeypatch):
    paths = concurrent_runs(tmp_path, 3)
    original = subprocess.run
    # Sequential execution would break the barrier at its timeout.
    started = threading.Barrier(3, timeout=10)

    def invoke(arguments, **kwargs):
        if len(arguments) == 2:
            started.wait()
        return original(arguments, **kwargs)

    monkeypatch.setattr(subprocess, "run", invoke)
    assert execute("hybridization", paths, concurrency=3) == [
        str(path.with_name(f"result{n}.h5")) for n, path in enumerate(paths, 1)]
    assert sorted(executed(commands)) == [["hybridization", str(path)] for path in paths]
    validations = [i for i, command in enumerate(commands) if "--validate" in command]
    assert max(validations) < commands.index(executed(commands)[0])


def test_concurrent_failure_waits_for_active_runs_and_starts_no_other(tmp_path, commands, monkeypatch):
    paths = concurrent_runs(tmp_path, 4)
    original = subprocess.run
    started = threading.Barrier(2, timeout=10)
    failing, finished = [], []

    def invoke(arguments, **kwargs):
        if len(arguments) == 2 and arguments[-1] == str(paths[0]):
            failing.append(threading.current_thread())
            started.wait()
            commands.append(list(arguments))
            raise subprocess.CalledProcessError(9, arguments)
        if len(arguments) == 2 and arguments[-1] == str(paths[1]):
            started.wait()
            # The failed run's worker exits only after recording its failure.
            failing[0].join(timeout=10)
            assert not failing[0].is_alive()
            finished.append(arguments[-1])
        return original(arguments, **kwargs)

    monkeypatch.setattr(subprocess, "run", invoke)
    with pytest.raises(subprocess.CalledProcessError) as failure:
        execute("hybridization", paths, concurrency=2)
    assert failure.value.returncode == 9
    assert finished == [str(paths[1])]
    assert sorted(executed(commands)) == [["hybridization", str(path)] for path in paths[:2]]


def test_concurrent_mpi_runs_accept_launcher_arguments(tmp_path, commands):
    paths = concurrent_runs(tmp_path, 2)
    execute("hybridization", paths, mpi=2, mpirun=["mpi launcher", "--bind-to", "none"], concurrency=2)
    assert sorted(executed(commands)) == [
        ["mpi launcher", "--bind-to", "none", "-np", "2", "hybridization", str(path)] for path in paths]


@pytest.mark.parametrize("arguments", [dict(concurrency=0), dict(concurrency=True),
                                       dict(mpi=1, mpirun=[])])
def test_invalid_concurrency_or_launcher_is_rejected_before_any_command(tmp_path, commands, arguments):
    paths = concurrent_runs(tmp_path, 2)
    with pytest.raises(ValueError):
        execute("hybridization", paths, **arguments)
    assert commands == []


def test_result_collisions_are_checked_for_existing_run_files(tmp_path, commands):
    paths = [tmp_path / "one.toml", tmp_path / "two.toml"]
    paths[0].write_text('[parameters]\ncount=1\n[output]\nresults="same.h5"\n')
    paths[1].write_text('[parameters]\ncount=2\n[output]\nresults="./same.h5"\n')
    with pytest.raises(ValueError, match="multiple runs"):
        execute("hybridization", paths)
    assert executed(commands) == []
    assert not any("--validate" in command for command in commands)


def test_existing_jobs_protect_inputs_of_every_run(tmp_path, commands):
    scientific = tmp_path / "data.h5"
    scientific.write_bytes(b"original measurements")
    first = write_run_file(tmp_path / "one.toml", SCHEMA,
                           parameters={"count": 1}, output={"results": "data.h5"})
    second = write_run_file(tmp_path / "two.toml", SCHEMA,
                            parameters={"count": 2}, input={"data": "data.h5"},
                            output={"results": "other.h5"})
    with pytest.raises(ValueError, match="overwrite input.data"):
        execute("hybridization", [first, second])
    assert executed(commands) == []
    assert not any("--validate" in command for command in commands)
    assert scientific.read_bytes() == b"original measurements"


def test_existing_jobs_check_additional_output_paths(tmp_path, commands, monkeypatch):
    schema = SCHEMA + '\n[output.final_omega]\ntype="path"\n'
    original = subprocess.run

    def invoke(arguments, **kwargs):
        result = original(arguments, **kwargs)
        if "--schema" in arguments:
            result.stdout = schema
        return result

    monkeypatch.setattr(subprocess, "run", invoke)
    paths = [write_run_file(tmp_path / f"run{n}.toml", schema,
                            parameters={"count": n},
                            output={"results": f"result{n}.h5", "final_omega": "shared.dat"})
             for n in (1, 2)]
    with pytest.raises(ValueError, match="multiple runs"):
        execute("dmft", paths)
    assert executed(commands) == []
    assert not any("--validate" in command for command in commands)


def test_each_runs_schema_protects_cross_run_sidecar_inputs(tmp_path, commands, monkeypatch):
    producer_schema = (SCHEMA + '\n[output.text]\ntype="bool"\ndefault=false\n'
                       '[output.text_directory]\ntype="path"\ndefault="."\n')
    consumer_schema = SCHEMA + '\n[input.initial_tau]\ntype="path"\n'
    sidecars = tmp_path / "sidecars"
    sidecars.mkdir()
    scientific = sidecars / "G_tau"
    scientific.write_bytes(b"original measurements")
    first = write_run_file(tmp_path / "producer.toml", producer_schema,
                           parameters={"count": 1},
                           output={"results": "first.h5", "text": True,
                                   "text_directory": "sidecars"})
    second = write_run_file(tmp_path / "consumer.toml", consumer_schema,
                            parameters={"count": 2}, input={"initial_tau": "sidecars/G_tau"},
                            output={"results": "second.h5"})
    original = subprocess.run

    def invoke(arguments, **kwargs):
        result = original(arguments, **kwargs)
        if "--schema" in arguments:
            result.stdout = producer_schema if arguments[-1] == str(first) else consumer_schema
        return result

    monkeypatch.setattr(subprocess, "run", invoke)
    with pytest.raises(ValueError, match="inside another run's output.text_directory"):
        execute("dmft", [first, second])
    assert [command[-1] for command in commands if "--schema" in command] == [str(first), str(second)]
    assert executed(commands) == []
    assert not any("--validate" in command for command in commands)
    assert scientific.read_bytes() == b"original measurements"


def test_existing_scheduler_xml_path_is_preserved(monkeypatch):
    recorded = []
    monkeypatch.setattr(tools, "check_existence", lambda _: None)
    monkeypatch.setattr(tools, "executeCommand", lambda args: recorded.append(args) or 0)
    assert tools.runApplication("sparsediag", "ed.in.xml", MPI=2, T=10, Tmax=30) == (0, "ed.out.xml")
    assert recorded == [["mpirun", "-np", "2", "sparsediag", "--mpi", "--Nmax", "1", "ed.in.xml",
                         "-T", "10", "--Tmax", "30"]]


@pytest.mark.parametrize("application, run", [("loop", "parm.in.xml"), ("/opt/alps/bin/spinmc", "parm.in.xml"),
                                              ("custom", "run.toml")])
def test_scheduler_launcher_refuses_toml_run_applications(monkeypatch, application, run):
    monkeypatch.setattr(tools, "executeCommand", lambda args: pytest.fail("must not launch"))
    with pytest.raises(ValueError, match="use pyalps.run_io.execute"):
        tools.runApplication(application, run)


def test_native_cthyb_launcher_roundtrip_when_cli_available(tmp_path, monkeypatch):
    executable = Path(alps_program("hybridization"))
    from pyalps import cthyb
    directory = tmp_path.resolve() / "scientific data with spaces"
    directory.mkdir()
    (directory / "delta.dat").write_text("0 -0.5 -0.5\n1 -0.5 -0.5\n2 -0.5 -0.5\n")
    parameters = {"BETA": 2., "N_ORBITALS": 2, "N_TAU": 2, "N_MEAS": 1,
                  "THERMALIZATION": 0, "SWEEPS": 10, "U": 1.}
    manifest = write_run_files(directory / "batch", [
        {"parameters": parameters, "input": {"delta": "delta.dat"},
         "output": {"results": f"result{n}.h5"}, "execution": {"seed": n}}
        for n in (1, 2)
    ], cthyb.schema())
    monkeypatch.chdir(tmp_path)
    results = execute(executable, manifest)
    assert results == [str(directory / "result1.h5"), str(directory / "result2.h5")]
    for filename in results:
        with hdf5.archive(filename, "r") as ar:
            assert ar["/run_config/application"] == "cthyb"
            assert ar.is_group("/simulation/results/Sign")


def test_native_hirschfye_launcher_and_modern_analysis_when_cli_available(tmp_path, monkeypatch):
    import numpy as np
    import pyalps
    from pyalps import alea

    executable = Path(alps_program("hirschfye"))
    directory = tmp_path / "Hirsch-Fye data with spaces"
    directory.mkdir()
    with hdf5.archive(directory / "g0.h5", "w") as archive:
        for flavor in range(2):
            archive[f"/G0_{flavor}"] = -2j / ((2*np.arange(4)+1)*np.pi)
    schema = subprocess.run([executable, "--schema"], check=True,
                            capture_output=True, text=True).stdout
    filename = directory / "result.h5"
    document = write_run_file(directory / "run.toml", schema,
        parameters={"BETA": 2., "U": 0., "N": 4, "NMATSUBARA": 4,
                    "SWEEPS": 37, "THERMALIZATION": 2,
                    "EPS_0": 0., "EPS_1": 0., "EPSSQ_0": 0., "EPSSQ_1": 0.},
        input={"g0": "g0.h5"}, output={"results": "result.h5"},
        execution={"bins": 8, "seed": 19})
    before = (directory / "g0.h5").read_bytes(), Path(document).read_bytes()
    monkeypatch.chdir(tmp_path)
    assert execute(executable, [document]) == [str(filename)]
    assert before == ((directory / "g0.h5").read_bytes(), Path(document).read_bytes())
    measured = {entry.props["observable"]: entry for entry in
                pyalps.loadMeasurements([str(filename)])[0]}
    assert set(measured) == {"Sign", "G_meas_up", "G_meas_down"}
    with hdf5.archive(filename) as archive:
        assert archive["/run_config/application"] == "hirschfye"
        for name, entry in measured.items():
            path = "/simulation/results/" + name
            assert archive[path + "/@kind"] == 5 and archive[path + "/@version"] == 1
            result = alea.BatchResult.read(archive, path)
            assert result.count == 37  # Warm-up includes the transition sweep.
            np.testing.assert_array_equal([value.mean for value in entry.y], result.mean)
            np.testing.assert_array_equal([value.error for value in entry.y], result.error)
            np.testing.assert_allclose(result.mean, 1. if name == "Sign" else -.5, atol=1e-13)
            np.testing.assert_allclose(result.error, 0., atol=1e-13)
            if name != "Sign":
                flavor = 0 if name == "G_meas_up" else 1
                np.testing.assert_array_equal(archive[f"/G_tau/{flavor}/mean/value"], result.mean)
                np.testing.assert_allclose(archive[f"/G_omega/{flavor}/mean/value"],
                                          -2j / ((2*np.arange(4)+1)*np.pi), atol=1e-13)
    assert {path.name for path in directory.iterdir()} == {"g0.h5", "run.toml", "result.h5"}


def test_native_dmft_ctint_python_workflow_and_u0_reference(tmp_path, monkeypatch):
    """An explicit Python job remains analysable after two real impurity runs."""
    dmft, ctint = Path(alps_program("dmft")), Path(alps_program("interaction"))
    import numpy as np
    from pyalps.load import loadDMFTIterations

    directory = tmp_path.resolve() / "scientific run with spaces"
    directory.mkdir()
    working = tmp_path.resolve() / "unrelated working directory"
    working.mkdir()
    # The driver finds its built-in solvers in the ALPS bin directory.
    monkeypatch.setenv("ALPS_BIN_PATH", str(ctint.parent))
    run_file = write_run_file(
        directory / "run.toml",
        parameters={"BETA": 2.0, "U": 0.0, "MU": 0.0, "H": 0.0, "t": 1.0, "N": 16,
                    "NMATSUBARA": 8, "FLAVORS": 2, "SITES": 1, "SWEEPS": 32,
                    "THERMALIZATION": 4, "CONVERGED": 0.0, "SYMMETRIZATION": True,
                    "ALPHA": -0.01, "MEASUREMENT_PERIOD": 1},
        output={"results": "results.h5"},
        execution={"solver": "interaction", "loop": "omega", "max_iterations": 2, "seed": 42})
    original_run = run_file.read_bytes()
    output = directory / "results.h5"
    monkeypatch.chdir(working)
    assert execute(dmft, run_file) == [str(output)]
    assert run_file.read_bytes() == original_run
    assert list(working.iterdir()) == []
    assert {path.name for path in directory.iterdir()} == {"run.toml", "results.h5"}

    # U=0 on the Bethe lattice has an exact fixed point. Check the physical
    # Matsubara convention, both flavors, and the complete iteration history.
    omega = (2 * np.arange(8) + 1) * np.pi / 2.0
    reference = -2j / (omega + np.sqrt(omega ** 2 + 4.0))
    with hdf5.archive(str(output), "r") as ar:
        saved = run_config.RunConfiguration()
        ar.set_context("/run_config")
        saved.load(ar)
        assert saved.application == "dmft"
        assert saved.execution["solver"] == "interaction"
        assert saved.output["results"] == str(output)
        assert saved.execution["seed"] == 42
        assert saved.execution["max_iterations"] == 2
        assert saved.parameters["BETA"] == 2.0 and saved.parameters["U"] == 0.0
        assert saved.parameters["SEMICIRCLE_HILBERT"] is False
        assert saved.origins["parameters.EPSSQ_0"] == "derived"
        assert saved.origins["parameters.EPS_0"] == "derived"
        assert "SEED" not in saved.parameters and "MAX_TIME" not in saved.parameters
        assert set(ar.list_children("/simulation/iteration")) == {"1", "2"}
        for iteration in (1, 2):
            base = f"/simulation/iteration/{iteration}/results"
            for flavor in (0, 1):
                frequency = np.asarray(ar[f"{base}/G_omega/{flavor}/mean/value"])
                bare = np.asarray(ar[f"{base}/G0_omega/{flavor}/mean/value"])
                time = np.asarray(ar[f"{base}/G_tau/{flavor}/mean/value"])
                assert frequency.shape == bare.shape == (8,)
                assert time.shape == (17,)
                assert np.iscomplexobj(frequency)
                assert np.all(np.isfinite(frequency)) and np.all(np.isfinite(time))
                np.testing.assert_allclose(frequency, reference, rtol=0, atol=1e-10)
                np.testing.assert_allclose(frequency, bare, rtol=0, atol=1e-12)
                np.testing.assert_allclose(time[0] + time[-1], -1.0, rtol=0, atol=1e-12)
                assert 0.0 <= -time[-1] <= 1.0
        for flavor in (0, 1):
            np.testing.assert_allclose(ar[f"/simulation/results/G_omega/{flavor}/mean/value"],
                                       reference, rtol=0, atol=1e-10)

    iterations = loadDMFTIterations([str(output)], observable="G_omega", measurements=["0", "1"])
    assert len(iterations) == 1 and len(iterations[0]) == 2
    assert {data[0].props["iteration"] for data in iterations[0]} == {"1", "2"}
    for flavors in iterations[0]:
        assert len(flavors) == 2
        assert {dataset.props["observable"] for dataset in flavors} == {"0", "1"}
        for dataset in flavors:
            assert dataset.props["filename"] == str(output)
            assert dataset.props["BETA"] == 2.0 and dataset.props["U"] == 0.0
            np.testing.assert_array_equal(dataset.x, np.arange(8))
            np.testing.assert_allclose(dataset.y, reference, rtol=0, atol=1e-10)
