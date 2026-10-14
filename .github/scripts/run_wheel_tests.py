"""Test the installed wheel and retain results across cibuildwheel interpreters."""

from pathlib import Path
import os
import platform
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[2]

# Regular manylinux checks retain installed SDK consumers but leave these source
# builds to native CI. Full compatibility/release runs exercise every contract
# in each distribution toolchain, including musl and macOS.
SOURCE_SDK_TESTS = (
    "test_embedded_defaults_and_mpi_isolation",
    "test_components_link_without_building_the_core_runtime",
    "test_utilities_without_platform_stacktrace",
    "test_make_intro_builds_without_sdk",
)


def main():
    import pyalps

    full_sdk = os.environ.get("PYALPS_TEST_FULL_SDK", "1")
    if full_sdk not in {"0", "1"}:
        raise ValueError("PYALPS_TEST_FULL_SDK must be 0 or 1")
    package = Path(pyalps.__file__).resolve()
    if package.is_relative_to(ROOT / "python"):
        raise RuntimeError(f"Wheel tests imported source-tree pyalps: {package}")
    reports = ROOT / "_build/wheel-reports"
    reports.mkdir(parents=True, exist_ok=True)
    identifier = f"{platform.python_implementation()}-{platform.python_version()}-{platform.machine()}"
    (reports / f"{identifier}-import.txt").write_text(f"Python: {sys.executable}\npyalps: {package}\n")
    # The C++ SDK is identical across wheel ABIs. Compile its consumers once
    # per platform; Python/downstream-extension tests still run on every ABI.
    suites = [str(ROOT / 'python/pyalps/tests')]
    selection = []
    if sys.version_info[:2] == (3, 12):
        if full_sdk == "1":
            suites.append(str(ROOT / 'tests/cmake'))
        else:
            sdk_tests = ROOT / 'tests/cmake/test_sdk_contracts.py'
            suites.extend((str(sdk_tests), str(ROOT / 'tests/cmake/test_python_sdk_version.py')))
            selection = [f"--deselect=tests/cmake/test_sdk_contracts.py::{name}"
                         for name in SOURCE_SDK_TESTS]
    elif sys.version_info[:2] > (3, 12):
        # Reused abi3 artifact: check interpreter-facing behavior, without
        # repeating its numerical regressions on every newer interpreter.
        suites = [str(ROOT / 'python/pyalps/tests' / name) for name in (
            'test_binding_surface.py', 'test_mapping_lifetimes.py', 'test_wheel_payload.py',
        )]
    with tempfile.TemporaryDirectory(prefix="alps-wheel-tests-") as work:
        try:
            return subprocess.call([
                sys.executable, "-m", "pytest", "-v", "--durations=20",
                f"--rootdir={ROOT}", "-o", "faulthandler_timeout=300", *suites, *selection,
                f"--junitxml={reports / f'{identifier}.xml'}",
            ], cwd=work)
        finally:
            # Report from inside the build container, including failed tests;
            # telemetry must never replace pytest's result.
            try:
                subprocess.call([
                    sys.executable, str(ROOT / ".github/scripts/report_ccache.py"),
                    "--phase", f"wheels-{identifier}",
                ], cwd=work)
            except OSError as error:
                print(f"Could not report compiler cache statistics: {error}", file=sys.stderr)


if __name__ == "__main__":
    raise SystemExit(main())
