# Copyright (c) 2026 Tobias Wolf
# SPDX-License-Identifier: MIT
#
# Adapted from skilledwolf/alpsx (private repository), branch
# agent/v3-execution-and-api-audit at d37e51227,
# python/qdk-chemistry-alps/tests/test_examples.py. Changes: the classic_dmrg
# marker replaces native_dmrg; the MPS and Q# examples are not carried over.
"""Run every example script.

The examples are the source for the code in the docs, so if one breaks, the
documentation is wrong. Each asserts its own physics internally; here we
only require that it exits cleanly and prints something.

They are sized to finish in seconds, so this runs with the rest of the
suite rather than as a nightly job.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

EXAMPLES = Path(__file__).resolve().parents[1] / "examples"

# Capabilities each example needs beyond the base stack (pyalps, qdk-chemistry).
REQUIREMENTS = {
    "quickstart.py": ["classic_dmrg"],
    "observables.py": ["classic_dmrg"],
    "large_active_space.py": ["classic_dmrg", "local_qdk"],
    "converging_dmrg.py": ["classic_dmrg"],
    "hardware_export.py": [],
    "interop.py": [],
    "molecule.py": ["pyscf"],
    "hubbard_reference.py": ["classic_dmrg"],
    "custom_spin_export.py": [],
    "fermion_pairing.py": [],
}


def _cases():
    for path in sorted(EXAMPLES.glob("*.py")):
        marks = [getattr(pytest.mark, name) for name in REQUIREMENTS[path.name]]
        yield pytest.param(path, marks=marks, id=path.stem)


def test_every_example_is_registered():
    """A new example must declare what it needs, or it silently never runs."""
    on_disk = {path.name for path in EXAMPLES.glob("*.py")}
    assert on_disk == set(REQUIREMENTS), (
        f"examples/ and REQUIREMENTS disagree: {on_disk ^ set(REQUIREMENTS)}"
    )


@pytest.mark.parametrize("script", list(_cases()))
def test_example_runs(script: Path):
    proc = subprocess.run(
        [sys.executable, str(script)],
        capture_output=True,
        text=True,
        timeout=900,
    )
    assert proc.returncode == 0, (
        f"{script.name} failed\n--- stdout ---\n{proc.stdout}\n"
        f"--- stderr ---\n{proc.stderr}"
    )
    assert proc.stdout.strip(), f"{script.name} printed nothing"
