"""Developer setup must isolate Python and prohibit native dependency builds."""

import importlib.util
import os
from pathlib import Path
import sys
import subprocess

import pytest


ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("alps_dev", ROOT / ".github/scripts/dev.py")
dev = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(dev)


def test_run_uses_the_managed_python_instead_of_the_callers_python(tmp_path):
    environment = tmp_path / "isolated"
    dev.venv.EnvBuilder().create(environment)
    scripts = environment / ("Scripts" if sys.platform == "win32" else "bin")
    result = dev.run(["python", "-c", "import sys; print(sys.prefix)"], capture=True,
                     env={**os.environ, "PATH": str(scripts) + os.pathsep + os.environ["PATH"]})
    assert Path(result.stdout.strip()).resolve() == environment.resolve()


@pytest.mark.parametrize("platform,expected", [("win-amd64", "x64"), ("win-arm64", "arm64")])
def test_windows_dependencies_match_python_not_host_cpu(monkeypatch, platform, expected):
    monkeypatch.setattr(dev.sysconfig, "get_platform", lambda: platform)
    assert dev.windows_architecture() == expected


def test_unsupported_python_architecture_has_clear_error(monkeypatch):
    monkeypatch.setattr(dev.sysconfig, "get_platform", lambda: "win32")
    with pytest.raises(RuntimeError, match="64-bit"):
        dev.windows_architecture()


def test_sdk_and_binding_configuration_disable_vcpkg_builds(tmp_path):
    sdk, bindings = dev.configuration(tmp_path / "dev", tmp_path / "dependencies", "arm64")
    assert "-DVCPKG_MANIFEST_MODE=OFF" in sdk
    assert "-DVCPKG_MANIFEST_MODE=OFF" in bindings
    assert "-DVCPKG_TARGET_TRIPLET=arm64-windows" in bindings
    assert any(option.startswith("-DCMAKE_TOOLCHAIN_FILE=") for option in bindings)
    assert "-DALPS_BUILD_TESTING=ON" in sdk and "-DALPS_ENABLE_MPI=OFF" in sdk


def test_unix_build_and_bindings_use_the_same_binary_environment(tmp_path):
    sdk, bindings = dev.configuration(tmp_path / "dev", tmp_path / "pixi")
    assert f"-DCMAKE_PREFIX_PATH={tmp_path / 'pixi'}" in sdk
    assert all(option in sdk for option in bindings)
    assert f"-DHDF5_ROOT={tmp_path / 'pixi'}" in bindings
    assert "-DBLA_VENDOR=Generic" in bindings


@pytest.mark.skipif(sys.platform == "win32" or not os.environ.get("PIXI_PROJECT_ROOT"),
                    reason="requires the managed Unix binary dependencies")
def test_managed_dependencies_ignore_a_competing_hdf5_config(tmp_path):
    prefix = Path(os.environ["CONDA_PREFIX"]).resolve()
    competitor = tmp_path / "system-hdf5"
    competitor.mkdir()
    (competitor / "hdf5-config.cmake").write_text(
        'message(FATAL_ERROR "Selected a competing system HDF5 config")\n')
    source = tmp_path / "source"
    source.mkdir()
    (source / "CMakeLists.txt").write_text('''
cmake_minimum_required(VERSION 3.27)
project(managed_dependencies LANGUAGES C CXX)
find_package(Boost 1.76 CONFIG REQUIRED COMPONENTS program_options)
find_package(HDF5 MODULE REQUIRED COMPONENTS C)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
get_target_property(boost_library Boost::program_options IMPORTED_LOCATION_RELEASE)
file(WRITE "${CMAKE_BINARY_DIR}/libraries.txt"
  "${boost_library}\\n${HDF5_C_LIBRARY_hdf5}\\n${BLAS_LIBRARIES}\\n${LAPACK_LIBRARIES}\\n")
''')
    _, shared = dev.configuration(tmp_path / "dev", prefix)
    build = tmp_path / "build"
    result = subprocess.run(["cmake", "-S", str(source), "-B", str(build), "-G", "Ninja",
                             f"-DHDF5_DIR={competitor}", *shared],
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    libraries = (build / "libraries.txt").read_text().splitlines()
    assert len(libraries) == 4
    for group in libraries:
        assert group
        for library in group.split(";"):
            assert Path(library).is_relative_to(prefix), library


def test_bootstrap_does_not_install_into_callers_python(tmp_path, monkeypatch):
    commands = []
    monkeypatch.setattr(dev, "run", lambda command, **kwargs: commands.append(command) or type("Result", (), {"returncode": 0})())
    monkeypatch.setattr(dev.venv.EnvBuilder, "create", lambda self, path: Path(path).mkdir(parents=True))
    monkeypatch.setattr(sys, "argv", ["dev.py", "configure"])
    assert dev.windows_environment(tmp_path, "x64") == 0
    assert str(commands[0][0]) == str(tmp_path / "venv/Scripts/python.exe")
    assert "--only-binary=:all:" in commands[0]
    assert str(commands[1][0]) == str(tmp_path / "venv/Scripts/python.exe")


def test_requirements_are_not_reinstalled_on_every_developer_command(tmp_path, monkeypatch):
    environment = tmp_path / "venv"
    (environment / "Scripts").mkdir(parents=True)
    (environment / "Scripts/python.exe").touch()
    (environment / "alps-requirements.sha256").write_text(dev.hashlib.sha256((ROOT / ".github/scripts/dev-requirements.txt").read_bytes()).hexdigest())
    commands = []
    monkeypatch.setattr(dev, "run", lambda command, **kwargs: commands.append(command) or type("Result", (), {"returncode": 0})())
    assert dev.windows_environment(tmp_path, "x64") == 0
    assert len(commands) == 1 and "pip" not in commands[0]


def test_x64_python_can_activate_a_native_arm64_host_compiler(tmp_path, monkeypatch):
    installer = tmp_path / "Microsoft Visual Studio/Installer/vswhere.exe"
    installer.parent.mkdir(parents=True)
    installer.touch()
    installation = tmp_path / "Visual Studio"
    vcvars = installation / "VC/Auxiliary/Build/vcvarsall.bat"
    vcvars.parent.mkdir(parents=True)
    vcvars.touch()
    monkeypatch.setenv("ProgramFiles(x86)", str(tmp_path))
    monkeypatch.delenv("VSCMD_ARG_TGT_ARCH", raising=False)
    monkeypatch.setattr(dev, "windows_host_architecture", lambda: "arm64")
    def execute(command, **kwargs):
        output = str(installation) if "-property" in command else "Path=native-compiler\nVSCMD_ARG_TGT_ARCH=x64\n"
        return subprocess.CompletedProcess(command, 0, stdout=output)
    monkeypatch.setattr(dev, "run", execute)
    environment = dev.msvc_environment("x64", tmp_path)
    assert 'arm64_x64' in (tmp_path / "msvc-env.cmd").read_text()
    assert environment["PATH"] == "native-compiler"


def test_dependency_update_refreshes_sdk_and_bindings_independently(tmp_path, monkeypatch):
    monkeypatch.setattr(dev, "ROOT", tmp_path)
    monkeypatch.setattr(dev.sys, "platform", "linux")
    monkeypatch.setenv("PIXI_PROJECT_ROOT", str(tmp_path))
    monkeypatch.setenv("PIXI_PLATFORM", "test-platform")
    monkeypatch.setenv("CONDA_PREFIX", str(tmp_path / "environment"))
    lock = tmp_path / "pixi.lock"
    lock.write_text("first dependencies")
    helper = tmp_path / ".github/scripts/dev.py"
    helper.parent.mkdir(parents=True)
    helper.write_text("first helper")
    commands = []
    def execute(command, **kwargs):
        commands.append([str(value) for value in command])
        if "-B" in command:
            Path(command[command.index("-B") + 1]).mkdir(parents=True, exist_ok=True)
        for option in command:
            if str(option).startswith("--config-settings=build-dir="):
                Path(str(option).split("=", 2)[2]).mkdir(parents=True, exist_ok=True)
        return subprocess.CompletedProcess(command, 0)
    monkeypatch.setattr(dev, "run", execute)
    def build(*options):
        commands.clear()
        monkeypatch.setattr(dev.sys, "argv", ["dev.py", *options])
        assert dev.main() == 0
        configure = next(command for command in commands if "-B" in command)
        bindings = [command for command in commands if "pip" in command]
        return "--fresh" in configure, bool(bindings and "--config-settings=cmake.args=--fresh" in bindings[0])
    assert build() == (True, True)
    assert build() == (False, False)
    lock.write_text("updated dependencies")
    assert build("--cpp-only") == (True, False)
    assert build() == (False, True)
    assert build() == (False, False)
    helper.write_text("updated dependency selection")
    assert build() == (True, True)
    assert build() == (False, False)
