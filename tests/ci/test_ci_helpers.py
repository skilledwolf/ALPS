"""Validate test reporting and installed-wheel execution."""

import importlib.util
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

import pytest


ROOT = Path(__file__).resolve().parents[2]


def load_helper(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / f".github/scripts/{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize("source,code", [
    ("def test_result():\n    assert True\n", 0),
    ("def test_result():\n    assert False\n", 1),
    ("# no tests\n", 5),
])
@pytest.mark.parametrize("version", [(3, 11), (3, 12), (3, 14)])
def test_wheel_runner_preserves_test_failures_and_empty_collection(tmp_path, monkeypatch, source, code, version):
    import types
    root = tmp_path / "checkout"
    for suite in ("python/pyalps/tests", "tests/cmake"):
        (root / suite).mkdir(parents=True)
    (root / "python/pyalps/tests/test_binding_surface.py").write_text(source)
    (root / "python/pyalps/tests/test_mapping_lifetimes.py").touch()
    (root / "python/pyalps/tests/test_wheel_payload.py").touch()
    scripts = root / ".github/scripts"
    scripts.mkdir(parents=True)
    # A failed telemetry helper must not mask either pytest failures or success.
    (scripts / "report_ccache.py").write_text("raise SystemExit(7)\n")
    monkeypatch.delenv("PYALPS_TEST_FULL_SDK", raising=False)
    helper = load_helper("run_wheel_tests")
    monkeypatch.setattr(helper, "sys", types.SimpleNamespace(version_info=version, executable=sys.executable))
    monkeypatch.setattr(helper, "ROOT", root)
    monkeypatch.setitem(sys.modules, "pyalps", types.SimpleNamespace(
        __file__=str(tmp_path / "installed/pyalps/__init__.py")))
    assert helper.main() == code
    reports = list((root / "_build/wheel-reports").glob("*.xml"))
    assert len(reports) == 1
    assert "test_binding_surface" in reports[0].read_text() or code == 5


@pytest.mark.parametrize("full_sdk", [None, "0", "1"])
@pytest.mark.parametrize("version", [(3, 11), (3, 12), (3, 14)])
def test_wheel_runner_keeps_installed_consumers_and_full_toolchain_coverage(
        tmp_path, monkeypatch, full_sdk, version):
    import types

    root = tmp_path / "checkout"
    python_tests = root / "python/pyalps/tests"
    sdk_tests = root / "tests/cmake"
    python_tests.mkdir(parents=True)
    sdk_tests.mkdir(parents=True)
    # The real package has its own pytest.ini. Explicit rootdir keeps SDK
    # deselection stable when pytest discovers that nested configuration.
    (python_tests.parent / "pytest.ini").write_text("[pytest]\n")
    (python_tests / "test_binding_surface.py").write_text(
        "def test_wheel_downstream_consumer(): pass\n"
        # An identically named wheel test must not be excluded accidentally.
        "def test_components_link_without_building_the_core_runtime(): pass\n")
    (python_tests / "test_mapping_lifetimes.py").touch()
    (python_tests / "test_wheel_payload.py").touch()
    (sdk_tests / "test_sdk_contracts.py").write_text(
        "import pytest\n"
        "def test_installed_sdk_preserves_parent_settings(): pass\n"
        "def test_relocated_sdk(): pass\n"
        "def test_embedded_defaults_and_mpi_isolation(): pass\n"
        "def test_components_link_without_building_the_core_runtime(): pass\n"
        "@pytest.mark.parametrize('mode', ['static', 'shared'])\n"
        "def test_utilities_without_platform_stacktrace(mode): pass\n"
        "def test_make_intro_builds_without_sdk(): pass\n")
    (sdk_tests / "test_python_sdk_version.py").write_text(
        "def test_python_rejects_mismatched_sdk(): pass\n")
    (sdk_tests / "test_build_configuration.py").write_text(
        "def test_embedded_build_keeps_parent_defaults(): pass\n")
    if full_sdk is None:
        monkeypatch.delenv("PYALPS_TEST_FULL_SDK", raising=False)
    else:
        monkeypatch.setenv("PYALPS_TEST_FULL_SDK", full_sdk)
    helper = load_helper("run_wheel_tests")
    monkeypatch.setattr(helper, "sys", types.SimpleNamespace(version_info=version, executable=sys.executable))
    monkeypatch.setattr(helper, "ROOT", root)
    monkeypatch.setitem(sys.modules, "pyalps", types.SimpleNamespace(
        __file__=str(tmp_path / "installed/pyalps/__init__.py")))
    assert helper.main() == 0
    report, = (root / "_build/wheel-reports").glob("*.xml")
    cases = ET.parse(report).findall(".//testcase")
    wheel_names = {case.attrib["name"] for case in cases
                   if "test_binding_surface" in case.attrib["classname"]}
    sdk_names = {case.attrib["name"] for case in cases
                 if "tests.cmake" in case.attrib["classname"]}
    assert wheel_names == {"test_wheel_downstream_consumer",
                           "test_components_link_without_building_the_core_runtime"}
    if version != (3, 12):
        assert not sdk_names
    else:
        assert {"test_installed_sdk_preserves_parent_settings", "test_relocated_sdk",
                "test_python_rejects_mismatched_sdk"} <= sdk_names
        source_names = {
            "test_embedded_defaults_and_mpi_isolation",
            "test_components_link_without_building_the_core_runtime",
            "test_utilities_without_platform_stacktrace[static]",
            "test_utilities_without_platform_stacktrace[shared]",
            "test_make_intro_builds_without_sdk",
            "test_embedded_build_keeps_parent_defaults",
        }
        if full_sdk == "0":
            assert not source_names & sdk_names
        else:
            assert source_names <= sdk_names


def test_wheel_runner_rejects_source_tree_import(monkeypatch):
    import types
    helper = load_helper("run_wheel_tests")
    monkeypatch.setitem(sys.modules, "pyalps", types.SimpleNamespace(
        __file__=str(ROOT / "python/pyalps/src/pyalps/__init__.py")))
    with pytest.raises(RuntimeError, match="source-tree"):
        helper.main()
