"""Test the installed wheel and retain results across cibuildwheel interpreters."""

from pathlib import Path
import platform
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[2]


def main():
    import pyalps

    package = Path(pyalps.__file__).resolve()
    if package.is_relative_to(ROOT / "python"):
        raise RuntimeError(f"Wheel tests imported source-tree pyalps: {package}")
    reports = ROOT / "_build/wheel-reports"
    reports.mkdir(parents=True, exist_ok=True)
    identifier = f"{platform.python_implementation()}-{platform.python_version()}-{platform.machine()}"
    (reports / f"{identifier}-import.txt").write_text(f"Python: {sys.executable}\npyalps: {package}\n")
    # Keep execution outside the checkout and pass explicit suite paths. The
    # wheel's SDK contract tests remain part of every selected interpreter.
    with tempfile.TemporaryDirectory(prefix="alps-wheel-tests-") as work:
        return subprocess.call([
            sys.executable, "-m", "pytest", "-v", "-o", "faulthandler_timeout=300",
            str(ROOT / "tests/pyalps"), str(ROOT / "tests/cmake"),
            f"--junitxml={reports / f'{identifier}.xml'}",
        ], cwd=work)


if __name__ == "__main__":
    raise SystemExit(main())
