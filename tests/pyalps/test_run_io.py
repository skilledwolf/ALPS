# Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
"""Python generation must round-trip through the application's native loader."""

import copy
from pathlib import Path

import numpy as np
import pytest

from pyalps import run_config
from pyalps.run_io import read_job_manifest, write_run_file, write_run_files


SCHEMA = '''
application = "writer-test"
schema_version = 1
[parameters.SWEEPS]
type = "int64"
required = true
min = 1
[parameters.J]
type = "float64"
default = 1.0
[parameters."MEASURE_AVERAGE[spin \\"quoted\\"]"]
type = "string"
[parameters.flags]
type = "bool[]"
[parameters.coefficients]
type = "complex128[]"
[parameters.offset]
type = "complex128"
[parameters.label]
type = "string"
[input.data]
type = "path"
[output.results]
type = "path"
required = true
[output.checkpoint]
type = "path"
[execution.seed]
type = "int64"
min = 0
max = 4294967295
'''


def test_round_trip_types_quoted_keys_paths_and_defaults(tmp_path):
    key = 'MEASURE_AVERAGE[spin "quoted"]'
    values = {
        "SWEEPS": np.int64(12), "flags": np.array([True, False]),
        "coefficients": np.array([1 + 2j, -3j]), "offset": 4 + 0j,
        "label": 'unicode α 🎲, slash \\, newline\n and "quotes"\x01\x7f', key: "Sz",
    }
    file = write_run_file(tmp_path / "nested" / "run.toml", SCHEMA,
                          parameters=values, input={"data": Path("data.h5")},
                          output={"results": "results.h5"})
    text = file.read_text(encoding="utf-8")
    assert "schema_version" not in text
    assert "application" not in text
    assert str(file.parent) not in text
    loaded = run_config.load(str(file), SCHEMA)
    assert loaded.parameters["SWEEPS"] == 12
    assert type(loaded.parameters["SWEEPS"]) is int
    assert loaded.parameters["J"] == 1.0
    assert loaded.origins["parameters.J"] == "default"
    assert loaded.parameters[key] == "Sz"
    assert loaded.parameters["label"] == values["label"]
    assert loaded.parameters["offset"] == 4 + 0j
    np.testing.assert_array_equal(loaded.parameters["flags"], values["flags"])
    np.testing.assert_array_equal(loaded.parameters["coefficients"], values["coefficients"])
    assert Path(loaded.input["data"]) == file.parent / "data.h5"
    assert Path(loaded.output["results"]) == file.parent / "results.h5"


def test_explicit_runs_vectors_and_nonmutating_seeds(tmp_path):
    runs = [
        {"parameters": {"SWEEPS": 2, "coefficients": [1j, 2j]},
         "output": {"results": "first.h5"}},
        {"parameters": {"SWEEPS": 3}, "execution": {"seed": 91},
         "output": {"results": "second.h5"}},
    ]
    original = copy.deepcopy(runs)
    manifest = write_run_files(tmp_path / 'batch "quoted"', runs, SCHEMA, baseseed=7)
    assert runs == original
    text = manifest.read_text(encoding="utf-8")
    assert text.count("[[runs]]") == 2
    assert 'application = "writer-test"' in text
    assert 'batch \\"quoted\\".task1.toml' in text
    first = run_config.load(str(tmp_path / 'batch "quoted".task1.toml'), SCHEMA)
    second = run_config.load(str(tmp_path / 'batch "quoted".task2.toml'), SCHEMA)
    np.testing.assert_array_equal(first.parameters["coefficients"], [1j, 2j])
    assert first.execution["seed"] == 7
    assert second.execution["seed"] == 91


def test_invalid_later_run_creates_no_files(tmp_path):
    runs = [{"parameters": {"SWEEPS": 2}, "output": {"results": "a.h5"}},
            {"parameters": {"SWEEPS": 0}, "output": {"results": "b.h5"}}]
    with pytest.raises((ValueError, RuntimeError)):
        write_run_files(tmp_path / "batch", runs, SCHEMA)
    assert list(tmp_path.iterdir()) == []


def test_refuse_overwrite_before_publishing_other_runs(tmp_path):
    existing = tmp_path / "batch.task2.toml"
    existing.write_text("owner's original data", encoding="utf-8")
    runs = [{"parameters": {"SWEEPS": 2}, "output": {"results": "a.h5"}},
            {"parameters": {"SWEEPS": 2}, "output": {"results": "b.h5"}}]
    with pytest.raises(FileExistsError):
        write_run_files(tmp_path / "batch", runs, SCHEMA)
    assert existing.read_text(encoding="utf-8") == "owner's original data"
    assert sorted(path.name for path in tmp_path.iterdir()) == [existing.name]


def test_explicit_overwrite_and_schema_errors(tmp_path):
    path = tmp_path / "run.toml"
    write_run_file(path, SCHEMA, parameters={"SWEEPS": 2}, output={"results": "a.h5"})
    with pytest.raises(FileExistsError):
        write_run_file(path, SCHEMA, parameters={"SWEEPS": 3}, output={"results": "b.h5"})
    write_run_file(path, SCHEMA, parameters={"SWEEPS": 3}, output={"results": "b.h5"},
                   overwrite=True)
    assert run_config.load(str(path), SCHEMA).parameters["SWEEPS"] == 3
    for invalid in ({"SWEEPS": "3"}, {"SWEEPS": 2, "unexpected": 1}):
        with pytest.raises((ValueError, RuntimeError, TypeError)):
            write_run_file(tmp_path / "invalid.toml", SCHEMA, parameters=invalid,
                           output={"results": "a.h5"})
    assert not (tmp_path / "invalid.toml").exists()


def test_seed_bounds_and_unknown_sections_are_not_silently_accepted(tmp_path):
    run = {"parameters": {"SWEEPS": 2}, "output": {"results": "a.h5"}}
    with pytest.raises((ValueError, RuntimeError)):
        write_run_files(tmp_path / "batch", [run, run], SCHEMA, baseseed=4294967295)
    with pytest.raises(ValueError, match="unknown run sections"):
        write_run_files(tmp_path / "batch", [dict(run, arbitrary={})], SCHEMA)
    with pytest.raises(TypeError, match="baseseed"):
        write_run_files(tmp_path / "batch", [run], SCHEMA, baseseed=True)
    with pytest.raises(ValueError, match="at least one"):
        write_run_files(tmp_path / "batch", [], SCHEMA)
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize("second_output", [
    {"results": "./same.h5"},
    {"results": "other.h5", "checkpoint": "same.h5"},
])
def test_reject_colliding_outputs_before_writing(tmp_path, second_output):
    runs = [{"parameters": {"SWEEPS": 2}, "output": {"results": "same.h5"}},
            {"parameters": {"SWEEPS": 3}, "output": second_output}]
    with pytest.raises(ValueError, match="multiple runs"):
        write_run_files(tmp_path / "batch", runs, SCHEMA)
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize("target", ["batch.task1.toml", "batch.task2.toml", "batch.job.toml"])
def test_outputs_cannot_overwrite_generated_configuration(tmp_path, target):
    runs = [{"parameters": {"SWEEPS": 2}, "output": {"results": target}},
            {"parameters": {"SWEEPS": 3}, "output": {"results": "other.h5"}}]
    with pytest.raises(ValueError, match="overwrite a run or job file"):
        write_run_files(tmp_path / "batch", runs, SCHEMA)
    assert list(tmp_path.iterdir()) == []


def test_reject_symlink_output_aliases(tmp_path):
    alias = tmp_path / "alias"
    try:
        alias.symlink_to(tmp_path, target_is_directory=True)
    except OSError:
        pytest.skip("directory symlinks unavailable")
    runs = [{"parameters": {"SWEEPS": 2}, "output": {"results": "same.h5"}},
            {"parameters": {"SWEEPS": 3}, "output": {"results": "alias/same.h5"}}]
    with pytest.raises(ValueError, match="multiple runs"):
        write_run_files(tmp_path / "batch", runs, SCHEMA)
    assert list(tmp_path.iterdir()) == [alias]


def test_schema_is_optional_and_omits_the_manifest_application(tmp_path):
    runs = [{"parameters": {"SWEEPS": 2}, "output": {"results": "a.h5"}},
            {"parameters": {"SWEEPS": 3}, "output": {"results": "b.h5"}}]
    manifest = write_run_files(tmp_path / "batch", runs)
    assert "application" not in manifest.read_text(encoding="utf-8")
    application, files = read_job_manifest(manifest)
    assert application is None
    assert [run_config.load(file, SCHEMA).parameters["SWEEPS"] for file in files] == [2, 3]
    with pytest.raises(ValueError, match="multiple runs"):
        write_run_files(tmp_path / "other", [runs[0], runs[0]])


def test_relative_outputs_are_anchored_at_the_run_file(tmp_path, monkeypatch):
    # data.h5 in the writer's directory is a different file from the run's data.h5.
    monkeypatch.chdir(tmp_path)
    path = write_run_file(tmp_path / "runs" / "run.toml", SCHEMA, parameters={"SWEEPS": 2},
                          input={"data": str(tmp_path / "data.h5")}, output={"results": "data.h5"})
    assert Path(run_config.load(path, SCHEMA).output["results"]) == path.parent / "data.h5"


def test_job_manifest_lists_runs_relative_to_the_manifest(tmp_path):
    manifest = tmp_path / "job.job.toml"
    manifest.write_text('application = "writer-test"\n[[runs]]\nfile = "run with spaces.toml"\n'
                        '[[runs]]\nfile = "second.toml"\n', encoding="utf-8")
    assert read_job_manifest(manifest) == (
        "writer-test", [(tmp_path / "run with spaces.toml").resolve(), (tmp_path / "second.toml").resolve()])


@pytest.mark.parametrize("contents, message", [
    ('application="a"\nruns=[]\n', "nonempty"),
    ('application=""\n[[runs]]\nfile="run.toml"\n', "application"),
    ('unknown=true\n[[runs]]\nfile="run.toml"\n', "unknown job-manifest"),
    ('[[runs]]\nfile="run.toml"\nextra=true\n', "only file"),
    ('[[runs]]\nfile="run.toml"\n[[runs]]\nfile="./run.toml"\n', "duplicate"),
])
def test_invalid_job_manifests_are_rejected(tmp_path, contents, message):
    manifest = tmp_path / "bad.job.toml"
    manifest.write_text(contents, encoding="utf-8")
    with pytest.raises(ValueError, match=message):
        read_job_manifest(manifest)
