# Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
# Use the canonical native example to keep the statistical workflow identical.
from pathlib import Path
import runpy

runpy.run_path(str(Path(__file__).resolve().parents[2] / "03-mc/01b-equilibration-and-convergence/tutorial1a.py"), run_name="__main__")
