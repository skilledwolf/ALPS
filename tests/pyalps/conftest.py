# Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
"""Locate the programs and build trees used by integration tests.

A program comes from ALPS_<NAME>_EXECUTABLE or the SDK named by ALPS_DIR
(<prefix>/share/alps, programs in <prefix>/bin). PATH and the programs bundled
in a wheel are not searched, so a test never picks up an unrelated ALPS
installation. ALPS_TEST_PROGRAMS=0 skips program tests while ALPS_DIR still
serves downstream builds; wheel tests use this instead of repeating the
source jobs' simulations for every interpreter.
Tutorial programs come from ALPS_TUTORIALS_BUILD_DIR, where each tutorial
project is built at its path relative to tutorials/.
A missing requirement skips its tests locally. CI sets
ALPS_REQUIRE_INTEGRATION=1 so that a misconfigured job fails instead of
passing by skipping.
"""
import os
from pathlib import Path
import shutil

import pytest


def unavailable(reason):
    if os.environ.get("ALPS_REQUIRE_INTEGRATION") == "1":
        pytest.fail(reason + " (ALPS_REQUIRE_INTEGRATION=1)", pytrace=False)
    pytest.skip(reason)


def alps_program(name):
    """Return the absolute path of an ALPS program, or skip/fail the test."""
    if os.environ.get("ALPS_TEST_PROGRAMS") == "0":
        pytest.skip("program tests are disabled by ALPS_TEST_PROGRAMS=0")
    variable = "ALPS_" + name.upper() + "_EXECUTABLE"
    explicit = os.environ.get(variable)
    if explicit:
        return str(Path(explicit).resolve(strict=True))
    sdk = os.environ.get("ALPS_DIR")
    found = sdk and shutil.which(name, path=str(Path(sdk).resolve().parents[1] / "bin"))
    if not found:
        unavailable(f"{name} is unavailable; set {variable} or ALPS_DIR")
    return str(Path(found).resolve())


def launch(executable, *runs):
    """Run TOML run files in one application process.

    pyalps.run_io.execute starts a schema query, a validation and the run for
    each file. Its own tests cover that path; the others need only the run.
    """
    import subprocess
    result = subprocess.run([str(executable), *map(str, runs)], capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, result.stdout + result.stderr


def tutorials_build():
    """Return ALPS_TUTORIALS_BUILD_DIR, which mirrors the tutorials/ layout."""
    value = os.environ.get("ALPS_TUTORIALS_BUILD_DIR")
    if not value:
        unavailable("set ALPS_TUTORIALS_BUILD_DIR to the built tutorial projects")
    return Path(value).resolve(strict=True)


@pytest.fixture
def openmp_examples():
    cache = tutorials_build() / "00-examples/CMakeCache.txt"
    if "ALPS_ENABLE_OPENMP:BOOL=ON" not in cache.read_text().splitlines():
        unavailable("configure the examples with ALPS_ENABLE_OPENMP=ON to test threaded updates")
