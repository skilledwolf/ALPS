"""File and programmatic runs share validation and archived provenance."""
# Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
from pathlib import Path
import subprocess
import sys
import numpy as np
import pytest
from pyalps import run_config, hdf5

SCHEMA = '''
application = "configuration-contract"
schema_version = 3
[parameters.count]
type = "int64"
required = true
min = 1
[parameters.rate]
type = "float64"
default = 0.25
[parameters.z]
type = "complex128"
[input.data]
type = "path"
[output.results]
type = "path"
default = "results.h5"
[execution.seed]
type = "int64"
default = 42
min = 0
'''


def test_native_sections_are_mappings_without_importing_mc_bindings():
    script = r'''
import sys
from pyalps import run_config
run = run_config.resolve('application="mapping"\nschema_version=1\n[parameters.count]\ntype="int64"\n[input]\n[output]\n[execution]', parameters={"count": 3})
assert dict(run.parameters) == {"count": 3}
assert "pyalps.ngs" not in sys.modules
'''
    result = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr


def test_native_file_and_mapping_agree(tmp_path):
    filename = tmp_path / "run.toml"
    filename.write_text('[parameters]\ncount=9007199254740993\n'
                        'z={real=1.0,imag=-2.0}\n[input]\ndata="data.h5"\n')
    source = {"count": 9007199254740993, "z": 1-2j}
    mapped = run_config.resolve(SCHEMA, parameters=source, input={"data": "data.h5"},
                                base_directory=tmp_path)
    loaded = run_config.load(filename, SCHEMA)
    assert loaded.source_file == str(filename.resolve())
    assert mapped.source_file == ""
    assert mapped.application == loaded.application == "configuration-contract"
    assert mapped.schema_version == loaded.schema_version == 3
    for section in ("parameters", "input", "output", "execution"):
        assert dict(getattr(mapped, section)) == dict(getattr(loaded, section))
    assert mapped.origins == loaded.origins
    assert mapped.origins["parameters.rate"] == "default"
    assert mapped.parameters["count"] == 9007199254740993
    assert Path(mapped.input["data"]) == tmp_path / "data.h5"
    assert source == {"count": 9007199254740993, "z": 1-2j}
    relative = run_config.resolve(SCHEMA, parameters=source, input={"data": "data.h5"})
    assert relative.input["data"] == "data.h5"


@pytest.mark.parametrize("sections, message", [
    ({"parameters": {"count": "1"}}, "count"),
    ({"parameters": {"count": True}}, "count"),
    ({"parameters": {"count": 0}}, "minimum"),
    ({"parameters": {"count": 1, "typo": 4}}, "unknown key"),
    ({"parameters": {"count": 1}, "execution": {"seed": -1}}, "minimum"),
    ({"parameters": {"count": 1}, "output": {"data": "x"}}, "unknown key"),
])
def test_programmatic_sections_are_validated(sections, message):
    with pytest.raises((ValueError, RuntimeError), match=message):
        run_config.resolve(SCHEMA, **sections)


@pytest.mark.parametrize("name", ["cthyb", "ctint", "maxent_c"])
def test_programmatic_solvers_reject_relative_input_output_aliases(tmp_path, monkeypatch, name):
    solver = pytest.importorskip("pyalps." + name)
    monkeypatch.chdir(tmp_path)
    if name == "cthyb":
        parameters = {"BETA": 2., "U": 1., "N_ORBITALS": 2, "N_TAU": 2,
                      "N_MEAS": 1, "THERMALIZATION": 0, "SWEEPS": 10}
        scientific = tmp_path / "delta.dat"
        scientific.write_text("0 -0.5 -0.5\n1 -0.5 -0.5\n2 -0.5 -0.5\n")
        inputs = {"delta": scientific.name}
    elif name == "ctint":
        parameters = {"BETA": 2., "U": 1., "MU": 0., "ALPHA": -0.01,
                      "N": 4, "NMATSUBARA": 4, "THERMALIZATION": 0, "SWEEPS": 10}
        scientific = tmp_path / "g0.h5"
        with hdf5.archive(str(scientific), "w") as ar:
            for flavor in (0, 1):
                ar[f"/G0_{flavor}"] = -2j / ((2 * np.arange(4) + 1) * np.pi)
        inputs = {"g0": scientific.name}
    else:
        parameters = {"BETA": 2., "NFREQ": 20, "OMEGA_MAX": 4.}
        scientific = tmp_path / "data.h5"
        with hdf5.archive(str(scientific), "w") as ar:
            ar["/Data"] = [-0.5, -0.3, -0.3, -0.5]
            ar["/Error"] = [0.01] * 4
        inputs = {"data": scientific.name}
    original = scientific.read_bytes()
    output = {"results": "./" + scientific.name}
    with pytest.raises((ValueError, RuntimeError), match="replace input"):
        solver.prepare(parameters, input=inputs, output=output)
    # A caller can construct or reload a deferred run without prepare(); solve
    # must enforce the same collision protection before any output is opened.
    deferred = run_config.resolve(solver.schema(), parameters=parameters, input=inputs, output=output)
    with pytest.raises((ValueError, RuntimeError), match="replace input"):
        solver.solve(deferred)
    assert scientific.read_bytes() == original


def test_prepared_run_keeps_paths_after_archive_reload_and_cwd_change(tmp_path, monkeypatch):
    solver = pytest.importorskip("pyalps.cthyb")
    directory = tmp_path.resolve() / "preparation"
    directory.mkdir()
    monkeypatch.chdir(directory)
    scientific = directory / "delta.dat"
    scientific.write_text("0 -0.5 -0.5\n1 -0.5 -0.5\n2 -0.5 -0.5\n")
    original = scientific.read_bytes()
    parameters = {"BETA": 2., "U": 1., "N_ORBITALS": 2, "N_TAU": 2,
                  "N_MEAS": 1, "THERMALIZATION": 0, "SWEEPS": 10}
    run = solver.prepare(parameters, input={"delta": "delta.dat"},
                         output={"results": "result.h5"})
    checkpoint = directory / "run.h5"
    with hdf5.archive(str(checkpoint), "w") as ar:
        ar["/run_config"] = run
    restored = run_config.RunConfiguration()
    with hdf5.archive(str(checkpoint), "r") as ar:
        ar.set_context("/run_config")
        restored.load(ar)
    elsewhere = tmp_path.resolve() / "execution"
    elsewhere.mkdir()
    monkeypatch.chdir(elsewhere)
    solver.solve(restored)
    assert restored.input["delta"] == str(scientific)
    assert restored.output["results"] == str(directory / "result.h5")
    assert scientific.read_bytes() == original
    assert list(elsewhere.iterdir()) == []
    with hdf5.archive(str(directory / "result.h5"), "r") as ar:
        assert ar.is_group("/simulation/results/Sign")
        ar.set_context("/run_config")
        recorded = run_config.RunConfiguration()
        recorded.load(ar)
        assert recorded.input["delta"] == str(scientific)
        assert recorded.output["results"] == str(directory / "result.h5")


def test_archived_run_roundtrip_and_transactional_failure(tmp_path):
    run = run_config.resolve(SCHEMA, parameters={"count": 2}, execution={"seed": 91})
    archive_file = str(tmp_path / "run.h5")
    with hdf5.archive(archive_file, "w") as ar:
        ar["/run_config"] = run
    restored = run_config.resolve(SCHEMA, parameters={"count": 99})
    with hdf5.archive(archive_file, "r") as ar:
        ar.set_context("/run_config")
        restored.load(ar)
    assert restored.application == run.application
    assert restored.origins == run.origins
    assert dict(restored.parameters) == dict(run.parameters)
    assert restored.execution["seed"] == 91
    with hdf5.archive(archive_file, "a") as ar:
        ar["/run_config/execution/format"] = "unsupported"
    with hdf5.archive(archive_file, "r") as ar:
        ar.set_context("/run_config")
        with pytest.raises((ValueError, RuntimeError), match="format"):
            restored.load(ar)
    assert restored.parameters["count"] == 2
    assert restored.execution["seed"] == 91


def test_native_formatter_roundtrips_scientific_values_and_quoted_keys(tmp_path):
    # Smoke test of the binding; the formatter's edge cases are tested natively.
    schema = '''
application = "formatter-contract"
schema_version = 1
[parameters]
"MEASURE[spin \\"α\\"]\\u0001" = {type="string"}
message = {type="string"}
count = {type="int64"}
ints = {type="int64[]"}
reals = {type="float64[]"}
flags = {type="bool[]"}
names = {type="string[]"}
empty = {type="string[]"}
z = {type="complex128"}
zs = {type="complex128[]"}
[input.data]
type = "string"
[output.results]
type = "string"
[execution.seed]
type = "int64"
'''
    key = 'MEASURE[spin "α"]\x01'
    values = {key: "Sz", "message": 'Unicode 🎲; \\\n\t\r"\'\'\'\x00\x01\x7f',
              "count": 9007199254740993, "ints": [-2**63, 2**63-1],
              "reals": [-0.0, np.nextafter(1., 2.), np.finfo(float).tiny, np.finfo(float).max],
              "flags": [True, False], "names": ['a,b', 'c\n"d'], "empty": [],
              "z": 1-2j, "zs": [complex(-0.0, 2.0), 3-4j]}
    run = run_config.resolve(schema, parameters=values, input={"data": "input.h5"},
                             output={"results": "results.h5"}, execution={"seed": 17})
    filename = tmp_path / "run.toml"
    filename.write_text(run_config.format(run), encoding="utf-8")
    loaded = run_config.load(filename, schema)
    assert loaded.application == run.application
    for section in ("parameters", "input", "output", "execution"):
        expected, actual = getattr(run, section), getattr(loaded, section)
        assert set(expected) == set(actual)
        for name in expected:
            np.testing.assert_array_equal(actual[name], expected[name])
    assert loaded.parameters["count"] == 9007199254740993
    assert np.signbit(loaded.parameters["reals"][0])
    assert np.signbit(loaded.parameters["zs"][0].real)
    empty = run_config.format(run_config.RunConfiguration())
    assert all(f"[{section}]" in empty for section in ("parameters", "input", "output", "execution"))
