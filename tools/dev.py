"""Build a usable ALPS checkout with binary dependencies and editable Python."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import sysconfig
import venv


ROOT = Path(__file__).resolve().parents[1]


def run(command, *, env=None, capture=False):
    command = [str(argument) for argument in command]
    if env is not None and sys.platform == "win32":
        # CreateProcess searches the parent's PATH when launching an executable;
        # passing a different child environment alone does not select our tools.
        command[0] = shutil.which(command[0], path=env.get("PATH")) or command[0]
    print("+ " + subprocess.list2cmdline(command), flush=True)
    return subprocess.run(command, cwd=ROOT, env=env, check=True,
                          capture_output=capture, text=True)


def windows_architecture():
    # Match the interpreter, including x64 Python running on an ARM64 host.
    architecture = {"win-amd64": "x64", "win-arm64": "arm64"}.get(sysconfig.get_platform())
    if not architecture:
        raise RuntimeError("Use 64-bit CPython 3.12+ for Windows x64 or ARM64.")
    return architecture


def windows_host_architecture():
    # Environment variables can describe the emulated x64 process on ARM64.
    # Ask Windows for the native machine so MSVC itself runs without emulation.
    import ctypes
    process = ctypes.c_ushort()
    native = ctypes.c_ushort()
    if not ctypes.windll.kernel32.IsWow64Process2(
        ctypes.c_void_p(-1), ctypes.byref(process), ctypes.byref(native)
    ):
        raise ctypes.WinError()
    return "arm64" if native.value == 0xAA64 else "x64"


def msvc_environment(architecture, directory):
    environment = dict(os.environ)
    if environment.get("VSCMD_ARG_TGT_ARCH", "").lower() == architecture and shutil.which("cl"):
        return environment
    installer = Path(os.environ.get("ProgramFiles(x86)", "C:/Program Files (x86)")) / "Microsoft Visual Studio/Installer/vswhere.exe"
    if not installer.is_file():
        raise RuntimeError("Install Visual Studio's Desktop development with C++ workload, then rerun this command.")
    installation = run([installer, "-latest", "-products", "*", "-property", "installationPath"], capture=True).stdout.strip()
    vcvars = Path(installation) / "VC/Auxiliary/Build/vcvarsall.bat"
    if not vcvars.is_file():
        raise RuntimeError("No MSVC C++ toolchain found. Install the Desktop development with C++ workload.")
    host = windows_host_architecture()
    target = ("arm64" if architecture == "arm64" else "arm64_x64") if host == "arm64" else (
        "amd64_arm64" if architecture == "arm64" else "amd64")
    if any(character in str(vcvars) for character in '"\r\n%'):
        raise RuntimeError("Unsupported characters in the Visual Studio installation path.")
    # One batch file activates the compiler and returns its environment. No
    # permanent PATH changes, registry writes, or compiler installation.
    script = directory / "msvc-env.cmd"
    script.write_text(f'@echo off\ncall "{vcvars}" {target} >nul\nif errorlevel 1 exit /b 1\nset\n')
    result = run([os.environ.get("COMSPEC", "cmd.exe"), "/d", "/c", str(script)], capture=True)
    for line in result.stdout.splitlines():
        key, separator, value = line.partition("=")
        if separator and key:
            environment[key.upper()] = value
    # Windows environment variable names are case insensitive; Python's dicts
    # are not. Keep only the activated PATH spelling.
    environment = {key.upper(): value for key, value in environment.items()}
    return environment


def windows_environment(directory, architecture):
    if sys.version_info < (3, 12):
        raise RuntimeError("The Windows bootstrap requires CPython 3.12 or newer.")
    environment = directory / "venv"
    python = environment / "Scripts/python.exe"
    if Path(sys.prefix).resolve() != environment.resolve():
        if not python.is_file():
            venv.EnvBuilder(with_pip=True).create(environment)
        requirements = ROOT / "tools/dev-requirements.txt"
        digest = hashlib.sha256(requirements.read_bytes()).hexdigest()
        stamp = environment / "alps-requirements.sha256"
        if not stamp.is_file() or stamp.read_text() != digest:
            run([python, "-m", "pip", "install", "--no-compile", "--only-binary=:all:", "-r", requirements])
            stamp.write_text(digest)
        return run([python, __file__, *sys.argv[1:]]).returncode
    return None


def binary_dependencies(directory, architecture):
    sys.path.insert(0, str(ROOT / ".github/scripts"))
    from ci_dependencies import download
    manifest = json.loads((ROOT / ".github/dependencies.json").read_text())
    package = "windows-" + architecture
    # An updated pin gets a new directory, preserving old builds until CMake
    # switches to the new verified dependencies. Failed downloads are not reused.
    digest = manifest["sha256"][package]
    destination = ROOT / "_build/dependencies" / f"{package}-{digest[:16]}"
    marker = destination / ".alps-dependency-sha256"
    if not marker.is_file() or marker.read_text() != digest:
        download(package, destination, manifest["repository"], manifest)
        marker.write_text(digest)
    return destination


def configuration(directory, dependency_prefix=None, architecture=None, testing=True):
    options = [
        "-G", "Ninja", "-DCMAKE_BUILD_TYPE=Release",
        f"-DALPS_BUILD_TESTING={'ON' if testing else 'OFF'}", "-DALPS_BUILD_APPLICATIONS=ON", "-DALPS_ENABLE_MPI=OFF",
        f"-DCMAKE_INSTALL_PREFIX={directory / 'install'}",
        "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON", "-DCMAKE_INSTALL_MESSAGE=LAZY",
    ]
    shared = []
    if architecture:
        shared = [
            f"-DCMAKE_TOOLCHAIN_FILE={dependency_prefix / 'scripts/buildsystems/vcpkg.cmake'}",
            f"-DVCPKG_INSTALLED_DIR={dependency_prefix / 'installed'}",
            f"-DVCPKG_TARGET_TRIPLET={architecture}-windows",
            "-DVCPKG_MANIFEST_MODE=OFF", "-DX_VCPKG_APPLOCAL_DEPS_SERIALIZED=ON",
        ]
        options += ["-DCMAKE_C_COMPILER=cl", "-DCMAKE_CXX_COMPILER=cl"]
    else:
        shared = [f"-DCMAKE_PREFIX_PATH={dependency_prefix}"]
    return options + shared, shared


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--jobs", type=int, default=min(os.cpu_count() or 2, 2),
                        help="compiler processes (default: 2, or 1 on a single-core host)")
    parser.add_argument("--cpp-only", action="store_true", help="skip the Python editable install")
    parser.add_argument("--fresh", action="store_true", help="reset CMake configuration after changing toolchains")
    parser.add_argument("command", nargs="?", choices=("dev", "configure", "test", "run"), default="dev")
    parser.add_argument("arguments", nargs=argparse.REMAINDER, help="command to execute with 'run'")
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    if args.command != "run" and args.arguments:
        parser.error("Put setup options before the command, for example: --jobs 2 test")
    if args.arguments[:1] == ["--"]:
        args.arguments.pop(0)
    windows = sys.platform == "win32"
    architecture = windows_architecture() if windows else None
    if not windows:
        project = os.environ.get("PIXI_PROJECT_ROOT", os.environ.get("PIXI_WORKSPACE_ROOT", ""))
        if not project or Path(project).resolve() != ROOT:
            pixi = shutil.which("pixi")
            if not pixi:
                raise RuntimeError("Install Pixi once (https://pixi.sh), then run 'pixi run --locked dev'.")
            return run([pixi, "run", "--locked", "python", __file__, *sys.argv[1:]]).returncode
    platform_name = "windows-" + architecture if windows else os.environ.get("PIXI_PLATFORM", sys.platform + "-" + platform.machine())
    directory = ROOT / "_build/dev" / platform_name
    directory.mkdir(parents=True, exist_ok=True)
    if windows:
        result = windows_environment(directory, architecture)
        if result is not None:
            return result
        environment = msvc_environment(architecture, directory)
        dependencies = binary_dependencies(directory, architecture)
        environment["VCPKG_ROOT"] = str(dependencies)
    else:
        environment = dict(os.environ)
        dependencies = Path(environment["CONDA_PREFIX"])
    scripts = Path(sys.executable).parent
    environment["PATH"] = os.pathsep.join((str(scripts), str(directory / "install/bin"), environment.get("PATH", "")))
    environment["ALPS_DIR"] = str(directory / "install/share/alps")
    environment["CMAKE_GENERATOR"] = "Ninja"
    environment["CMAKE_BUILD_PARALLEL_LEVEL"] = str(args.jobs)
    environment["CTEST_PARALLEL_LEVEL"] = str(min(args.jobs, 2))
    if shutil.which("ccache", path=environment["PATH"]):
        environment["CMAKE_C_COMPILER_LAUNCHER"] = "ccache"
        environment["CMAKE_CXX_COMPILER_LAUNCHER"] = "ccache"
    configure, shared = configuration(directory, dependencies, architecture, testing=args.command == "test")
    environment["ALPS_TEST_CMAKE_ARGS"] = json.dumps(["-G", "Ninja", *shared])
    if args.command == "run":
        if not args.arguments:
            parser.error("run needs a command, for example: run python -c \"import pyalps\"")
        if not (directory / "install/share/alps/ALPSConfig.cmake").is_file():
            raise RuntimeError("Build this checkout first: python tools/dev.py")
        return run(args.arguments, env=environment).returncode
    # CMake otherwise preserves package locations from the previous dependency
    # set. SDK and bindings have separate stamps because either can build alone.
    inputs = ([ROOT / ".github/dependencies.json", ROOT / "tools/dev-requirements.txt"]
              if windows else [ROOT / "pixi.lock"])
    fingerprint = hashlib.sha256(b"\0".join(path.read_bytes() for path in inputs)).hexdigest()
    stamps = {name: directory / name / "alps-environment.sha256" for name in ("sdk", "bindings")}
    refresh = {name: args.fresh or not stamp.is_file() or stamp.read_text() != fingerprint
               for name, stamp in stamps.items()}
    if refresh["sdk"]:
        configure.append("--fresh")
    run(["cmake", "-S", ROOT, "-B", directory / "sdk", *configure], env=environment)
    stamps["sdk"].write_text(fingerprint)
    if args.command == "configure":
        return 0
    run(["cmake", "--build", directory / "sdk"], env=environment)
    run(["cmake", "--install", directory / "sdk"], env=environment)
    if not args.cpp_only:
        command = [sys.executable, "-m", "pip", "install", "--no-build-isolation", "--no-deps",
                   "-e", str(ROOT / "python/pyalps"),
                   f"--config-settings=build-dir={directory / 'bindings'}"]
        command += ["--config-settings=cmake.define." + option.removeprefix("-D") for option in shared]
        if refresh["bindings"]:
            command.append("--config-settings=cmake.args=--fresh")
        run(command, env=environment)
        stamps["bindings"].write_text(fingerprint)
        run([sys.executable, "-c", "import pyalps; print('Editable pyalps:', pyalps.__file__)"], env=environment)
    if args.command == "test":
        run(["ctest", "--test-dir", directory / "sdk", "--output-on-failure", "--no-tests=error"], env=environment)
        if not args.cpp_only:
            environment["PYALPS_TEST_DOWNSTREAM_EXPORT"] = "1"
            run([sys.executable, "-m", "pytest", "tests/pyalps", "tests/cmake", "-q"], env=environment)
    print(f"ALPS is ready. SDK: {directory / 'install'}\n"
          "Run programs with: python tools/dev.py run <command>\n"
          "After C++ changes, rerun this command; Python edits are already live.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (RuntimeError, OSError, subprocess.CalledProcessError) as error:
        raise SystemExit(f"ALPS development setup failed: {error}") from error
