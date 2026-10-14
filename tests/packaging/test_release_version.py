"""Regression coverage for release version drift and stale build artifacts."""

import importlib.util
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import zipfile

from packaging.version import Version
import pytest

SCRIPT = Path(__file__).resolve().parents[2] / ".github/scripts" / "check_release_version.py"
SPEC = importlib.util.spec_from_file_location("check_release_version", SCRIPT)
release = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(release)


@pytest.fixture
def versions(tmp_path):
    def write(core="3.0.0", python="3.0.0"):
        (tmp_path / "cmake").mkdir(exist_ok=True)
        (tmp_path / "cmake/ALPS_VERSION.txt").write_text(core + "\n")
        project_dir = tmp_path / "python/pyalps"
        project_dir.mkdir(parents=True, exist_ok=True)
        (project_dir / "pyproject.toml").write_text(
            f'[project]\nname = "pyalps"\nversion = "{python}"\n'
        )
        return tmp_path
    return write


@pytest.mark.parametrize("ref,expected", [
    ("refs/heads/master", "3.0.0"),
    ("refs/tags/v3.0.0", "3.0.0"),
    ("refs/tags/v3.0.0-beta.2", "3.0.0b2"),
    ("refs/tags/v3.0.0-rc.1", "3.0.0rc1"),
])
def test_standalone_dynamic_version(versions, monkeypatch, ref, expected):
    root = versions()
    monkeypatch.delenv("ALPS_VERSION_PRERELEASE", raising=False)
    (root / "python/pyalps/pyproject.toml").write_text(
        '[project]\nname = "pyalps"\ndynamic = ["version"]\n'
    )
    assert release.check_version(root, ref) == Version(expected)


def test_dynamic_version_rejects_stale_tag_and_conflicting_label(versions, monkeypatch):
    root = versions()
    (root / "python/pyalps/pyproject.toml").write_text(
        '[project]\nname = "pyalps"\ndynamic = ["version"]\n'
    )
    with pytest.raises(ValueError, match="disagrees with ALPS_VERSION.txt"):
        release.check_version(root, "refs/tags/v2.3.4")
    monkeypatch.setenv("ALPS_VERSION_PRERELEASE", "beta.1")
    with pytest.raises(ValueError, match="disagrees with ALPS_VERSION_PRERELEASE"):
        release.check_version(root, "refs/tags/v3.0.0")


def test_prerelease_sdist_keeps_its_version_without_the_build_environment(tmp_path):
    """Exercise the backend: a user rebuild must keep the published version."""
    repository = SCRIPT.parents[2]
    project = repository / "python/pyalps"
    core = (repository / "cmake/ALPS_VERSION.txt").read_text().strip()
    environment = {**os.environ, "GITHUB_REF": f"refs/tags/v{core}-beta.2"}
    environment.pop("ALPS_VERSION_PRERELEASE", None)
    subprocess.run(
        [sys.executable, "-c",
         "from scikit_build_core.build import build_sdist; "
         "import sys; build_sdist(sys.argv[1])", str(tmp_path)],
        cwd=project, env=environment, check=True, capture_output=True, text=True,
    )
    with tarfile.open(tmp_path / f"pyalps-{core}b2.tar.gz") as archive:
        prefix = f"pyalps-{core}b2/"
        required = {
            "ALPS_VERSION.txt", "LICENSE.txt", "CMakeLists.txt", "pyproject.toml",
            "_build_support/alps_version.py", "src/pyalps/__init__.py",
            "cpp/hdf5.cpp", "cpp/maxent_c.cpp", "cpp/cthyb.cpp", "cpp/ctint.cpp",
            "_build_support/runtime_manifest.py", "_build_support/pyalpsConfig.cmake",
        }
        assert {prefix + name for name in required} <= set(archive.getnames())
        archive.extractall(tmp_path, filter="data")
    environment.pop("GITHUB_REF")
    unpacked = tmp_path / f"pyalps-{core}b2"
    completed = subprocess.run(
        [sys.executable, "-c",
         "from scikit_build_core.build import prepare_metadata_for_build_wheel; "
         "print(prepare_metadata_for_build_wheel('metadata'))"],
        cwd=unpacked, env=environment, check=True, capture_output=True, text=True,
    )
    assert f"pyalps-{core}b2.dist-info" in completed.stdout


@pytest.mark.parametrize("ref", [
    "", "refs/heads/master", "refs/heads/v4.0.0", "refs/pull/142/merge", "refs/tags/v3.0.0"
])
def test_matching_versions(versions, ref):
    assert release.check_version(versions(), ref) == Version("3.0.0")


def test_ci_rejects_original_release_using_github_ref(versions):
    result = subprocess.run(
        [sys.executable, release.__file__, "--root", str(versions("2.3.4", "2.3.4b1"))],
        env={**os.environ, "GITHUB_REF": "refs/tags/v3.0.0"},
        capture_output=True,
        text=True,
    )
    assert result.returncode == 1
    assert "Release tag v3.0.0 disagrees" in result.stderr


@pytest.mark.parametrize("ref", [
    "", "refs/heads/master", "refs/pull/142/merge", "refs/tags/v3.0.0"
])
def test_sdk_and_python_must_agree_even_on_branches(versions, ref):
    with pytest.raises(ValueError, match="disagrees with ALPS_VERSION.txt"):
        release.check_version(versions("2.3.4", "3.0.0"), ref)


@pytest.mark.parametrize("label,suffix", [
    ("alpha.1", "a1"), ("beta.2", "b2"), ("rc.3", "rc3"), ("dev.4", ".dev4")
])
def test_prerelease_tags_match_pep440_versions(versions, label, suffix):
    python = "3.0.0" + suffix
    version = release.check_version(versions(python=python), "refs/tags/v3.0.0-" + label)
    assert version == Version(python)


@pytest.mark.parametrize("python,tag", [
    ("3.0.0b1", "v3.0.0"), ("3.0.0", "v3.0.0-beta.1"), ("3.0.0b1", "v3.0.0-beta.2")
])
def test_prerelease_cannot_be_published_as_final_or_different_prerelease(versions, python, tag):
    with pytest.raises(ValueError, match="disagrees with pyproject.toml"):
        release.check_version(versions(python=python), "refs/tags/" + tag)


@pytest.mark.parametrize("tag", [
    "v3.0", "v3.0.0-beta", "v3.0.0-final", "v3.0.0_1", "v03.0.0", "v3.0.0+local"
])
def test_malformed_tags_fail_closed(versions, tag):
    with pytest.raises(ValueError, match="Invalid release tag"):
        release.check_version(versions(), "refs/tags/" + tag)


@pytest.mark.parametrize("core", ["3.0", "v3.0.0", "3.0.0-beta.1", "3.0.0\n2.3.4", "03.0.0"])
def test_invalid_numeric_core(versions, core):
    with pytest.raises(ValueError, match="must contain MAJOR.MINOR.PATCH"):
        release.check_version(versions(core=core), "")


@pytest.fixture
def artifact(tmp_path):
    def write(
        kind, version="3.0.0", metadata_version=None, name="pyalps",
        tags="cp313-cp313-manylinux_2_28_x86_64", directory=None,
    ):
        directory = tmp_path if directory is None else directory
        directory.mkdir(parents=True, exist_ok=True)
        metadata = (
            f"Metadata-Version: 2.1\nName: {name}\nVersion: {metadata_version or version}\n"
        ).encode()
        if kind == "wheel":
            path = directory / f"{name}-{version}-{tags}.whl"
            with zipfile.ZipFile(path, "w") as archive:
                archive.writestr(f"{name}-{version}.dist-info/METADATA", metadata)
        else:
            path = directory / f"{name}-{version}.tar.gz"
            with tarfile.open(path, "w:gz") as archive:
                member = tarfile.TarInfo(f"{name}-{version}/PKG-INFO")
                member.size = len(metadata)
                archive.addfile(member, io.BytesIO(metadata))
        return path
    return write


@pytest.mark.parametrize("version", ["3.0.0", "3.0.0b1"])
def test_matching_wheel_and_sdist(tmp_path, artifact, version):
    artifact("wheel", version)
    artifact("sdist", version)
    release.check_distributions(tmp_path, Version(version))


def test_stale_artifact_rejects_the_entire_batch(tmp_path, artifact):
    artifact("wheel")
    artifact("sdist", "2.3.4b1")
    with pytest.raises(ValueError, match="expected a pyalps 3.0.0"):
        release.check_distributions(tmp_path, Version("3.0.0"))


@pytest.mark.parametrize("kind", ["wheel", "sdist"])
def test_renaming_an_old_artifact_does_not_fix_its_version(tmp_path, artifact, kind):
    artifact(kind, metadata_version="2.3.4b1")
    with pytest.raises(ValueError, match="metadata does not describe pyalps 3.0.0"):
        release.check_distributions(tmp_path, Version("3.0.0"))


def test_other_project_is_rejected(tmp_path, artifact):
    artifact("wheel", name="other")
    with pytest.raises(ValueError, match="expected a pyalps 3.0.0"):
        release.check_distributions(tmp_path, Version("3.0.0"))


def test_empty_dist_is_rejected(tmp_path):
    with pytest.raises(ValueError, match="No distributions"):
        release.check_distributions(tmp_path, Version("3.0.0"))


def test_unexpected_file_is_rejected(tmp_path):
    (tmp_path / "README.txt").write_text("not a distribution")
    with pytest.raises(ValueError, match="Unexpected distribution"):
        release.check_distributions(tmp_path, Version("3.0.0"))


@pytest.fixture
def release_batch(tmp_path, artifact):
    def write(profile="full"):
        paths = {"sdist": artifact("sdist", directory=tmp_path / "cibw-sdist")}
        platforms = {
            "manylinux": "manylinux_2_28_x86_64",
            "musllinux": "musllinux_1_2_x86_64",
            "macos": "macosx_15_0_arm64",
        }
        for family, platform in platforms.items():
            if profile == "minimal" and family != "manylinux":
                continue
            for interpreter, abi in (("cp310", "cp310"), ("cp311", "cp311"), ("cp312", "abi3")):
                paths[f"{family}/{interpreter}"] = artifact(
                    "wheel", tags=f"{interpreter}-{abi}-{platform}",
                    directory=tmp_path / f"cibw-wheels-{family}",
                )
        return paths
    return write


def matrix_policy(profile="full"):
    return release.load_matrix_policy(SCRIPT.parents[2], profile)


@pytest.mark.parametrize("profile,count", [("minimal", 4), ("full", 10)])
def test_complete_release_matrix(tmp_path, release_batch, profile, count):
    paths = release_batch(profile)
    verified = release.check_distributions(tmp_path, Version("3.0.0"), matrix_policy(profile))
    assert set(verified) == set(paths.values())
    assert len(verified) == count


@pytest.mark.parametrize("missing", ["manylinux/cp310", "musllinux/cp311", "macos/cp312"])
def test_missing_release_variant_is_rejected(tmp_path, release_batch, missing):
    release_batch()[missing].unlink()
    with pytest.raises(ValueError, match="Missing release wheels"):
        release.check_distributions(tmp_path, Version("3.0.0"), matrix_policy())


def test_missing_release_sdist_is_rejected(tmp_path, release_batch):
    release_batch()["sdist"].unlink()
    with pytest.raises(ValueError, match="Expected exactly one sdist, found 0"):
        release.check_distributions(tmp_path, Version("3.0.0"), matrix_policy())


def test_duplicate_filenames_are_checked_before_flattening(tmp_path, release_batch):
    original = release_batch()["manylinux/cp310"]
    duplicate = tmp_path / "cibw-unexpected"
    duplicate.mkdir()
    shutil.copyfile(original, duplicate / original.name)
    with pytest.raises(ValueError, match="Duplicate distribution filename"):
        release.check_distributions(tmp_path, Version("3.0.0"), matrix_policy())


def test_second_normalized_sdist_is_rejected(tmp_path, release_batch, artifact):
    release_batch()
    artifact("sdist", name="PyALPS", directory=tmp_path / "cibw-extra-sdist")
    with pytest.raises(ValueError, match="Expected exactly one sdist, found 2"):
        release.check_distributions(tmp_path, Version("3.0.0"), matrix_policy())


@pytest.mark.parametrize("tags", [
    "cp312-cp312-manylinux_2_28_x86_64",
    "cp313-abi3-manylinux_2_28_x86_64",
    "cp312-abi3-manylinux_2_28_aarch64",
    "py3-none-any",
    "cp312-abi3-manylinux_2_28_x86_64.linux_x86_64",
])
def test_unexpected_release_tags_are_rejected(tmp_path, release_batch, artifact, tags):
    release_batch()
    artifact("wheel", tags=tags, directory=tmp_path / "cibw-extra-wheel")
    with pytest.raises(ValueError, match="unexpected wheel tag"):
        release.check_distributions(tmp_path, Version("3.0.0"), matrix_policy())


def test_two_wheels_covering_one_variant_are_rejected(tmp_path, release_batch, artifact):
    release_batch()
    artifact("wheel", tags="1-cp312-abi3-manylinux_2_28_x86_64", directory=tmp_path / "cibw-extra-wheel")
    with pytest.raises(ValueError, match="Conflicting wheels"):
        release.check_distributions(tmp_path, Version("3.0.0"), matrix_policy())


def test_compound_repair_platform_tags_are_one_variant(tmp_path, release_batch, artifact):
    release_batch()["manylinux/cp312"].unlink()
    artifact(
        "wheel", tags="cp312-abi3-manylinux_2_17_x86_64.manylinux2014_x86_64",
        directory=tmp_path / "cibw-wheels-manylinux",
    )
    release.check_distributions(tmp_path, Version("3.0.0"), matrix_policy())


def test_wheel_cannot_cover_two_platform_families(tmp_path, release_batch, artifact):
    paths = release_batch()
    paths["manylinux/cp312"].unlink()
    paths["musllinux/cp312"].unlink()
    artifact(
        "wheel", tags="cp312-abi3-manylinux_2_28_x86_64.musllinux_1_2_x86_64",
        directory=tmp_path / "cibw-wheels-manylinux",
    )
    with pytest.raises(ValueError, match="overlapping release variants"):
        release.check_distributions(tmp_path, Version("3.0.0"), matrix_policy())


@pytest.mark.parametrize("profile,families", [
    ("minimal", ["manylinux"]),
    ("full", ["manylinux", "musllinux", "macos"]),
])
def test_matrix_cli_emits_only_json(profile, families):
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--matrix", profile],
        capture_output=True, text=True, check=True,
    )
    matrix = json.loads(result.stdout)
    assert [platform["family"] for platform in matrix] == families
    assert all(set(platform) == {"os", "target", "arch", "family"} for platform in matrix)


@pytest.mark.parametrize("complete", [True, False])
def test_release_cli_stages_only_complete_artifacts(tmp_path, versions, release_batch, complete):
    paths = release_batch()
    if not complete:
        paths["musllinux/cp311"].unlink()
    # Keep project metadata outside the recursive artifact tree.
    downloads = tmp_path / "downloads"
    downloads.mkdir()
    for path in list(tmp_path.glob("cibw-*")):
        shutil.move(path, downloads / path.name)
    root = versions()
    (root / ".github").mkdir()
    shutil.copyfile(SCRIPT.parents[2] / ".github/wheel-matrix.json", root / ".github/wheel-matrix.json")
    destination = tmp_path / "dist"
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--root", str(root), "--ref", "refs/tags/v3.0.0",
         "--dist", str(downloads), "--matrix-policy", "full", "--stage-dist", str(destination)],
        capture_output=True, text=True,
    )
    if complete:
        assert result.returncode == 0, result.stderr
        assert {path.name for path in destination.iterdir()} == {path.name for path in paths.values()}
    else:
        assert result.returncode == 1
        assert "Missing release wheels" in result.stderr
        assert not destination.exists()


def test_staging_refuses_existing_output(tmp_path, release_batch):
    artifacts = list(release_batch().values())
    destination = tmp_path / "dist"
    destination.mkdir()
    stale = destination / "stale.whl"
    stale.write_text("stale")
    with pytest.raises(FileExistsError):
        release.stage_distributions(artifacts, destination)
    assert list(destination.iterdir()) == [stale]
