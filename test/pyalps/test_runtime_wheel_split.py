#!/usr/bin/env python3
# Copyright (C) 2026 by the ALPS collaboration
# SPDX-License-Identifier: MIT

"""Tests for the release-time pyalps runtime-wheel splitter."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import zipfile

import pytest


REPOSITORY = Path(__file__).resolve().parents[2]
SPLITTER_PATH = REPOSITORY / ".github" / "scripts" / "split_pyalps_runtime_wheels.py"
SPEC = importlib.util.spec_from_file_location("split_pyalps_runtime_wheels", SPLITTER_PATH)
assert SPEC is not None and SPEC.loader is not None
splitter = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = splitter
SPEC.loader.exec_module(splitter)


def _fake_wheel(
    directory: Path,
    python_tag: str,
    *,
    platform_tag: str = "manylinux_2_28_x86_64",
    runtime_marker: bytes = b"same-runtime",
) -> Path:
    version = "2.3.4b1"
    dist_info = f"pyalps-{version}.dist-info"
    path = directory / f"pyalps-{version}-{python_tag}-{python_tag}-{platform_tag}.whl"
    members = {
        "pyalps/__init__.py": b"# Python-specific package\n",
        f"pyalps/_ext/example.{python_tag}.so": python_tag.encode(),
        "pyalps/bin/spinmc": runtime_marker + b"-program",
        "pyalps/lib/libalps.so.2": runtime_marker + b"-alps",
        "pyalps/xml/ALPS.xsl": runtime_marker + b"-xml",
        "pyalps.libs/libhdf5-deadbeef.so.1": runtime_marker + b"-hdf5",
        f"{dist_info}/METADATA": (
            "Metadata-Version: 2.1\n"
            "Name: pyalps\n"
            f"Version: {version}\n"
            "Requires-Dist: numpy>=1.26\n\n"
        ).encode(),
        f"{dist_info}/WHEEL": (
            "Wheel-Version: 1.0\n"
            "Root-Is-Purelib: false\n"
            f"Tag: {python_tag}-{python_tag}-{platform_tag}\n\n"
        ).encode(),
        f"{dist_info}/licenses/LICENSE.txt": b"MIT\n",
        f"{dist_info}/RECORD": b"input record is intentionally replaced\n",
    }
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as wheel:
        for name, data in members.items():
            wheel.writestr(name, data)
    return path


def test_split_creates_one_runtime_wheel_per_platform(tmp_path):
    input_dir = tmp_path / "input"
    output_dir = tmp_path / "output"
    input_dir.mkdir()
    _fake_wheel(input_dir, "cp310")
    _fake_wheel(input_dir, "cp314")

    outputs = splitter.split_wheel_directory(input_dir, output_dir)

    assert len(outputs) == 3
    runtime = output_dir / (
        "pyalps_runtime-2.3.4b1-py3-none-manylinux_2_28_x86_64.whl"
    )
    assert runtime in outputs
    with zipfile.ZipFile(runtime) as wheel:
        names = set(wheel.namelist())
        assert "pyalps/bin/spinmc" in names
        assert "pyalps/lib/libalps.so.2" in names
        assert "pyalps/xml/ALPS.xsl" in names
        assert "pyalps.libs/libhdf5-deadbeef.so.1" in names
        assert "pyalps/__init__.py" not in names
        wheel_metadata = wheel.read("pyalps_runtime-2.3.4b1.dist-info/WHEEL")
        assert b"Tag: py3-none-manylinux_2_28_x86_64" in wheel_metadata

    for python_tag in ("cp310", "cp314"):
        main = output_dir / (
            f"pyalps-2.3.4b1-{python_tag}-{python_tag}-manylinux_2_28_x86_64.whl"
        )
        with zipfile.ZipFile(main) as wheel:
            names = set(wheel.namelist())
            assert "pyalps/__init__.py" in names
            assert not any(splitter._is_runtime_member(name) for name in names)
            metadata = wheel.read("pyalps-2.3.4b1.dist-info/METADATA")
            assert b"Requires-Dist: pyalps-runtime==2.3.4b1" in metadata

    for output in outputs:
        splitter.verify_record(output)


def test_split_rejects_runtime_that_differs_between_abis(tmp_path):
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    _fake_wheel(input_dir, "cp310")
    _fake_wheel(input_dir, "cp314", runtime_marker=b"changed-runtime")

    with pytest.raises(ValueError, match="runtime payload differs"):
        splitter.split_wheel_directory(input_dir, tmp_path / "output")
