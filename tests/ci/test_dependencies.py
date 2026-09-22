"""Binary downloads must verify integrity before extracting anything."""

import hashlib
import importlib.util
import io
import json
from pathlib import Path
import tarfile
import re

import pytest


ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("ci_dependencies", ROOT / ".github/scripts/ci_dependencies.py")
dependencies = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(dependencies)


def archive(name="include/example.h"):
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w:gz") as stream:
        entry = tarfile.TarInfo(name)
        entry.size = 6
        stream.addfile(entry, io.BytesIO(b"header"))
    return buffer.getvalue()


def install(data, tmp_path, monkeypatch, checksum=None):
    monkeypatch.setattr(dependencies.urllib.request, "urlopen", lambda *args, **kwargs: io.BytesIO(data))
    manifest = {"release": "ci-dependencies-test", "sha256": {
        "test": checksum or hashlib.sha256(data).hexdigest()}}
    dependencies.download("test", tmp_path / "installed", "owner/repository", manifest)


def test_verified_archive_is_installed(tmp_path, monkeypatch):
    install(archive(), tmp_path, monkeypatch)
    assert (tmp_path / "installed/include/example.h").read_bytes() == b"header"


def test_corrupt_archive_is_never_installed(tmp_path, monkeypatch):
    with pytest.raises(ValueError, match="Checksum mismatch"):
        install(archive(), tmp_path, monkeypatch, "0" * 64)
    assert not (tmp_path / "installed").exists()


def test_archive_cannot_escape_destination(tmp_path, monkeypatch):
    with pytest.raises(tarfile.FilterError):
        install(archive("../escaped"), tmp_path, monkeypatch)
    assert not (tmp_path / "installed").exists()
    assert not (tmp_path / "escaped").exists()


def test_missing_pin_fails_without_network(tmp_path, monkeypatch):
    def unexpected_download(*args, **kwargs):
        pytest.fail("Missing pins must fail before any download")
    monkeypatch.setattr(dependencies.urllib.request, "urlopen", unexpected_download)
    with pytest.raises(ValueError, match="No pinned binary"):
        dependencies.download("missing", tmp_path / "installed", "owner/repository",
                              {"release": "ci-dependencies-test", "sha256": {}})


def test_pins_cover_every_source_and_packaging_platform():
    spec = importlib.util.spec_from_file_location("ci_matrix", ROOT / ".github/scripts/ci_matrix.py")
    matrix = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(matrix)
    source = json.loads((ROOT / ".github/ci-matrix.json").read_text())
    builds = matrix.select_matrix(source, "full")
    required = {build["dependency"] for build in builds["include"]}
    required.update(("windows-x64", "windows-arm64", f"boost-{source['boost_default']}-manylinux-x64"))
    manifest = json.loads((ROOT / ".github/dependencies.json").read_text())
    assert manifest["sha256"].keys() == required
    assert all(re.fullmatch(r"[a-f0-9]{64}", checksum) for checksum in manifest["sha256"].values())
