# Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
from pathlib import Path
import runpy

example = Path(__file__).resolve().parents[2] / "03-mc/01-autocorrelations/tutorial1a.py"
runpy.run_path(str(example))["run"](update="cluster", thermalization=1000,
                                  sweeps=100000, prefix="parm1b")
