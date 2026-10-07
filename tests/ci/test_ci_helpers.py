"""Validate conservative CI selection and retained test evidence."""

import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest


ROOT = Path(__file__).resolve().parents[2]


def load_helper(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / f".github/scripts/{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_source_matrix_expands_every_manifest_entry():
    manifest = json.loads((ROOT / ".github/ci-matrix.json").read_text())
    builds = load_helper("ci_matrix").select_matrix(manifest)["include"]
    assert len(builds) == len(manifest["builds"])
    for entry, build in zip(manifest["builds"], builds):
        assert entry.items() <= build.items()
        assert build["boost"] in json.loads((ROOT / ".github/dependencies.json").read_text())["boost"]


def test_matrix_cli_writes_github_outputs(tmp_path):
    output = tmp_path / "output"
    summary = tmp_path / "summary"
    subprocess.run(
        [sys.executable, str(ROOT / ".github/scripts/ci_matrix.py")],
        env={**os.environ, "GITHUB_OUTPUT": str(output),
             "GITHUB_STEP_SUMMARY": str(summary),
             "GITHUB_EVENT_NAME": "workflow_dispatch"},
        check=True, capture_output=True, text=True,
    )
    values = dict(line.split("=", 1) for line in output.read_text().splitlines())
    builds = json.loads(values["matrix"])["include"]
    manifest = json.loads((ROOT / ".github/ci-matrix.json").read_text())
    assert len(builds) == len(manifest["builds"])
    assert len({build["id"] for build in builds}) == len(builds)
    assert str(len(builds)) in summary.read_text()


@pytest.mark.parametrize("problem", ["duplicate", "checksum", "empty"])
def test_reject_broken_matrix(problem):
    manifest = json.loads((ROOT / ".github/ci-matrix.json").read_text())
    dependencies = json.loads((ROOT / ".github/dependencies.json").read_text())
    if problem == "duplicate":
        manifest["builds"].append(manifest["builds"][0])
    elif problem == "checksum":
        dependencies["boost"]["1.76.0"] = "not a checksum"
    else:
        manifest["builds"] = []
    with pytest.raises(ValueError):
        load_helper("ci_matrix").select_matrix(manifest, dependencies=dependencies)


def test_mixed_junit_results(tmp_path):
    report = tmp_path / "results.xml"
    report.write_text(
        '<testsuites><testsuite><testcase name="pass"/>'
        '<testcase name="failure"><failure>bad result</failure></testcase>'
        '<testcase name="error"><error>setup failed</error></testcase>'
        '<testcase name="skip"><skipped/></testcase></testsuite></testsuites>'
    )
    result = load_helper("junit_summary").summarize([report])
    assert "| 1 | 2 | 1 |" in result


def test_missing_junit_is_explicit():
    result = load_helper("junit_summary").summarize([])
    assert "No test reports" in result


@pytest.mark.parametrize("event", ["push", "schedule", "workflow_dispatch", None])
def test_non_pr_events_always_keep_complete_coverage(event):
    assert load_helper("ci_matrix").select_profile(event, ["README.md"]) == "full"


@pytest.mark.parametrize("paths,profile", [
    (["README.md", "docs/testing.md"], "fast"),
    (["tutorials/12-optical-lattice/01-bandstructure/bandstructure.py",
      "tests/tutorials/test_bandstructure.py"], "fast"),
    (["python/pyalps/src/pyalps/hdf5.py"], "representative"),
    (["src/alps/hdf5/archive.cpp", "tests/pyalps/test_archive_dtypes.py"], "representative"),
    (["src/alps/hdf5/archive.hpp"], "full"),
    (["src/alps/hdf5/CMakeLists.txt"], "full"),
    (["cmake/ALPSTesting.cmake"], "full"),
    ([".github/workflows/build.yml"], "full"),
    (["python/pyalps/pyproject.toml"], "full"),
    (["tests/conftest.py"], "full"),
    (["pytest.ini"], "full"),
    (["third_party/new-dependency/header.hpp"], "full"),
    (["new-component/code.cpp"], "full"),
    (["README.md", "CMakePresets.json"], "full"),
    ([], "full"),
    (None, "full"),
])
def test_pr_profile_is_conservative(paths, profile):
    assert load_helper("ci_matrix").select_profile("pull_request", paths) == profile


def test_representative_matrix_keeps_distinct_supported_contracts():
    manifest = json.loads((ROOT / ".github/ci-matrix.json").read_text())
    builds = load_helper("ci_matrix").select_matrix(manifest, "representative")["include"]
    assert len(builds) < len(manifest["builds"])
    assert all(build["pr"] for build in builds)
    assert any(build["cc"] == "gcc-11" for build in builds)
    assert any(build["cc"].startswith("clang-") for build in builds)
    assert any(build["os"].startswith("macos-") for build in builds)
    assert any(build["boost"] == "1.76.0" for build in builds)
    assert any(build["extras"] for build in builds)
    assert any(build["extensive"] for build in builds)
    assert all(build["mpi"] == "ON" for build in builds)


def test_complete_manifest_is_validated_even_for_reduced_prs():
    manifest = json.loads((ROOT / ".github/ci-matrix.json").read_text())
    unused = next(entry for entry in manifest["builds"] if not entry.get("pr"))
    manifest["builds"].append(unused)
    with pytest.raises(ValueError, match="duplicate"):
        load_helper("ci_matrix").select_matrix(manifest, "representative")


@pytest.mark.parametrize("profile", ["full", "representative", "fast"])
def test_packaging_keeps_all_families_and_supported_interpreter_boundaries(profile):
    wheels, smoke = load_helper("ci_matrix").packaging_matrices(profile)
    assert {entry["family"] for entry in wheels["plat"]} == {"manylinux", "musllinux", "macos"}
    assert smoke["python"][0] == "3.11" and smoke["python"][-1] == "3.14"
    if profile == "full":
        assert smoke["python"] == ["3.11", "3.12", "3.13", "3.14"]
        assert any(entry["os"] == "macos-26" for entry in smoke["plat"])
        assert all("{11,12,13,14}" in entry["build"] for entry in wheels["plat"])
    else:
        assert all("{11,14}" in entry["build"] for entry in wheels["plat"])


def test_missing_pr_diff_selects_complete_matrix(tmp_path):
    output = tmp_path / "output"
    subprocess.run(
        [sys.executable, str(ROOT / ".github/scripts/ci_matrix.py")],
        env={**os.environ, "GITHUB_OUTPUT": str(output),
             "GITHUB_EVENT_NAME": "pull_request", "GITHUB_EVENT_PATH": str(tmp_path / "missing"),
             "GITHUB_STEP_SUMMARY": str(tmp_path / "summary")},
        check=True, capture_output=True, text=True,
    )
    values = dict(line.split("=", 1) for line in output.read_text().splitlines())
    assert values["profile"] == "full"
    assert values["native"] == "true"


def test_changed_paths_include_both_sides_of_a_rename(tmp_path, monkeypatch):
    # Moving native code into the fast path must still exercise native coverage.
    def git(*args):
        return subprocess.check_output(["git", *args], cwd=tmp_path, text=True).strip()
    git("init", "-q")
    git("config", "user.name", "CI selector test")
    git("config", "user.email", "test@example.invalid")
    (tmp_path / "src").mkdir()
    (tmp_path / "src/native.cpp").write_text("native source\n")
    git("add", ".")
    git("commit", "-qm", "base")
    base = git("rev-parse", "HEAD")
    (tmp_path / "tests/tutorials").mkdir(parents=True)
    git("mv", "src/native.cpp", "tests/tutorials/moved.py")
    git("commit", "-qm", "move")
    head = git("rev-parse", "HEAD")
    event = tmp_path / "event.json"
    event.write_text(json.dumps({"pull_request": {"base": {"sha": base}, "head": {"sha": head}}}))
    helper = load_helper("ci_matrix")
    monkeypatch.setattr(helper, "ROOT", tmp_path)
    paths = helper.changed_paths({"GITHUB_EVENT_NAME": "pull_request", "GITHUB_EVENT_PATH": str(event)})
    assert set(paths) == {"src/native.cpp", "tests/tutorials/moved.py"}
    assert helper.select_profile("pull_request", paths) != "fast"


def test_junit_summary_includes_elapsed_time(tmp_path):
    report = tmp_path / "results.xml"
    report.write_text('<testsuite><testcase name="a" time="1.25"/>'
                      '<testcase name="b" time="0.5"/></testsuite>')
    result = load_helper("junit_summary").summarize([report])
    assert "Seconds" in result
    assert "| 2 | 0 | 0 | 1.75 |" in result


def test_standalone_directory_does_not_hide_new_native_build_files():
    helper = load_helper("ci_matrix")
    assert helper.select_profile("pull_request", ["tests/tutorials/CMakeLists.txt"]) == "full"
    assert helper.select_profile("pull_request", ["tutorials/12-optical-lattice/01-bandstructure/helper.hpp"]) == "full"


@pytest.mark.parametrize("source,code", [
    ("def test_result():\n    assert True\n", 0),
    ("def test_result():\n    assert False\n", 1),
    ("# no tests\n", 5),
])
def test_wheel_runner_preserves_test_failures_and_empty_collection(tmp_path, monkeypatch, source, code):
    import types
    root = tmp_path / "checkout"
    for suite in ("pyalps", "cmake"):
        (root / "tests" / suite).mkdir(parents=True)
    (root / "tests/pyalps/test_contract.py").write_text(source)
    helper = load_helper("run_wheel_tests")
    monkeypatch.setattr(helper, "ROOT", root)
    monkeypatch.setitem(sys.modules, "pyalps", types.SimpleNamespace(
        __file__=str(tmp_path / "installed/pyalps/__init__.py")))
    assert helper.main() == code
    reports = list((root / "_build/wheel-reports").glob("*.xml"))
    assert len(reports) == 1
    assert "test_contract" in reports[0].read_text() or code == 5


def test_wheel_runner_rejects_source_tree_import(monkeypatch):
    import types
    helper = load_helper("run_wheel_tests")
    monkeypatch.setitem(sys.modules, "pyalps", types.SimpleNamespace(
        __file__=str(ROOT / "python/pyalps/src/pyalps/__init__.py")))
    with pytest.raises(RuntimeError, match="source-tree"):
        helper.main()
