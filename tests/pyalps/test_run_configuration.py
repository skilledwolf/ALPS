"""File and programmatic runs share validation and archived provenance."""
# Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
from pathlib import Path
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


def test_native_file_and_mapping_agree(tmp_path):
    filename = tmp_path / "run.toml"
    filename.write_text('[parameters]\ncount=9007199254740993\n'
                        'z={real=1.0,imag=-2.0}\n[input]\ndata="data.h5"\n')
    source = {"count": 9007199254740993, "z": 1-2j}
    mapped = run_config.resolve(SCHEMA, parameters=source, input={"data": "data.h5"},
                                base_directory=tmp_path)
    loaded = run_config.load(filename, SCHEMA)
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
