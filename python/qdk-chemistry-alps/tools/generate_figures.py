# Copyright (c) 2026 Tobias Wolf
# SPDX-License-Identifier: MIT
#
# Adapted from skilledwolf/alpsx (private repository), branch
# agent/v3-execution-and-api-audit at d37e51227,
# python/qdk-chemistry-alps/tools/generate_figures.py. Changes: the runs use
# alps_classic_dmrg instead of alps_dmrg; the CI-memory figure uses that
# calculator's estimate (1 kB per determinant) instead of an MPS
# reconstruction workspace; the convergence figure is in units of t without
# the chemical-accuracy line; figures get an opaque white background instead
# of currentColor text, as plain Markdown loads them as images.
"""Regenerate the data figures used in the documentation.

    python python/qdk-chemistry-alps/tools/generate_figures.py

Each figure is written twice: an SVG for the docs and a JSON file holding
the plotted values. The figures are committed rather than built with the
docs, which keeps matplotlib and the solver stack out of reading them.
Needs matplotlib and the solver stack of the examples.
"""

from __future__ import annotations

import json
import logging
import os
from math import comb
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
# The font stack below is resolved by the viewer, not here; matplotlib only
# needs it for layout measurement and warns about every name it lacks.
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

FIGURES = Path(
    os.environ.get(
        "QDK_ALPS_FIGURE_DIR",
        Path(__file__).resolve().parents[1] / "docs" / "figures",
    )
)

INK = "#2a2a2a"
ACCENT = "#0e9384"      # teal
CONTRAST = "#c2632a"    # warm
MUTED = "#8a9a99"

plt.rcParams.update({
    "svg.fonttype": "none",           # keep text as text, not paths
    "svg.hashsalt": "qdk-chemistry-alps",  # stable element ids between runs
    "font.family": "system-ui, -apple-system, Segoe UI, Helvetica, Arial, sans-serif",
    "font.size": 9,
    "axes.edgecolor": INK,
    "axes.labelcolor": INK,
    "text.color": INK,
    "xtick.color": INK,
    "ytick.color": INK,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
})


def _save(fig, name: str, data: dict) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    svg = FIGURES / f"{name}.svg"
    fig.savefig(svg, format="svg", bbox_inches="tight", metadata={"Date": None})
    plt.close(fig)
    svg.write_text("\n".join(line.rstrip() for line in svg.read_text().splitlines()) + "\n")
    (FIGURES / f"{name}.json").write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    print(f"  wrote {svg.name} and {name}.json")


def _calculator(max_states: int):
    from qdk_chemistry.algorithms import create

    import qdk_chemistry_alps  # noqa: F401

    calc = create("multi_configuration_calculator", "alps_classic_dmrg")
    calc.settings().set("max_bond_dim", max_states)
    return calc


def _hubbard(n: int, U: float):
    from qdk_chemistry.data import LatticeGraph
    from qdk_chemistry.utils.model_hamiltonians import create_hubbard_hamiltonian

    return create_hubbard_hamiltonian(LatticeGraph.chain(n), epsilon=0.0, t=1.0, U=U)


# --------------------------------------------------------------------------
# 1. DMRG convergence -- mirrors examples/converging_dmrg.py
# --------------------------------------------------------------------------

def convergence() -> dict:
    states = [4, 8, 16, 32, 64, 128]
    hamiltonian = _hubbard(8, 8.0)
    energies = [float(_calculator(m).run(hamiltonian, 4, 4)[0]) for m in states]

    best = min(energies)
    errors = [e - best for e in energies]

    fig, ax = plt.subplots(figsize=(5.2, 3.2))
    # The best point has zero error and cannot be drawn on a log axis.
    ax.loglog(states[:-1], errors[:-1], "o-", color=ACCENT, linewidth=1.8,
              markersize=5, label="error above best")
    ax.scatter([states[-1]], [errors[-2] / 30], marker="v", color=ACCENT, s=40,
               zorder=3)
    ax.annotate("best", xy=(states[-1], errors[-2] / 30),
                xytext=(-6, 10), textcoords="offset points",
                fontsize=8, color=ACCENT, ha="right")
    ax.set_xlabel("kept states (max_bond_dim)")
    ax.set_ylabel("energy above best (units of t)")
    ax.set_xticks(states)
    ax.set_xticklabels([str(m) for m in states])
    ax.grid(True, which="major", color=MUTED, alpha=0.25, linewidth=0.6)

    _save(fig, "dmrg-convergence",
          {"max_bond_dim": states, "energies": energies, "model": "hubbard L=8 U=8 (4,4)",
           "calculator": "alps_classic_dmrg"})
    return {"states": states, "energies": energies}


# --------------------------------------------------------------------------
# 2. The CI-vector limit -- why the wavefunction representation has to change
# --------------------------------------------------------------------------

def ci_wall() -> dict:
    orbitals = list(range(4, 25, 2))
    determinants = [comb(n, n // 2) ** 2 for n in orbitals]
    # The calculator's estimate: 1 kB per determinant of the sector (about
    # 460 B measured for qdk-chemistry's objects, plus the engine's arrays).
    memory_gb = [d * 1e3 / 1e9 for d in determinants]
    guard = 4.0

    fig, ax = plt.subplots(figsize=(5.2, 3.2))
    ax.semilogy(orbitals, memory_gb, "o-", color=ACCENT, linewidth=1.8,
                markersize=4.5)
    ax.axhline(guard, color=CONTRAST, linestyle="--", linewidth=1.2)
    ax.text(4.2, guard * 1.5, "ci_memory_limit_gb (4 GB)", color=CONTRAST, fontsize=8)

    crossing = next(n for n, gb in zip(orbitals, memory_gb, strict=True) if gb > guard)
    ax.axvspan(crossing - 1, orbitals[-1] + 0.5, color=CONTRAST, alpha=0.07)
    ax.text((crossing + orbitals[-1]) / 2, 1e-4, "RDMs returned\ninstead",
            ha="center", fontsize=8, color=CONTRAST)
    ax.text((orbitals[0] + crossing) / 2, 1e-4, "CI vector\nreturned",
            ha="center", fontsize=8, color=ACCENT)

    ax.set_xlabel("active orbitals (half filled)")
    ax.set_ylabel("estimated CI-vector memory (GB)")
    ax.set_xticks(orbitals[::2])
    ax.set_xlim(orbitals[0] - 0.5, orbitals[-1] + 0.5)
    ax.grid(True, which="major", color=MUTED, alpha=0.25, linewidth=0.6)

    _save(fig, "ci-wall",
          {"orbitals": orbitals, "determinants": determinants, "memory_gb": memory_gb,
           "bytes_per_determinant": 1000, "guard_gb": guard})
    return {"crossing": crossing}


# --------------------------------------------------------------------------
# 3. Spin-spin correlations -- what the output actually looks like
# --------------------------------------------------------------------------

def correlations() -> dict:
    n = 8
    calc = _calculator(128)
    calc.run(_hubbard(n, 4.0), n // 2, n // 2)
    spin = np.asarray(calc.spin_spin_correlations()).real

    fig, (ax, bx) = plt.subplots(1, 2, figsize=(8.6, 3.1),
                                 gridspec_kw={"width_ratios": [1, 1.05],
                                              "wspace": 0.45})

    limit = np.max(np.abs(spin - np.diag(np.diag(spin))))
    image = ax.imshow(spin, cmap="RdBu_r", vmin=-limit, vmax=limit)
    ax.set_title(r"$\langle S_i \cdot S_j \rangle$", fontsize=9)
    ax.set_xlabel("orbital $j$")
    ax.set_ylabel("orbital $i$")
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    bar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    bar.outline.set_edgecolor(INK)
    bar.ax.tick_params(color=INK, labelcolor=INK)

    distance = np.arange(1, n)
    averaged = [float(np.mean([spin[i, i + d] for i in range(n - d)])) for d in distance]
    bx.axhline(0.0, color=MUTED, linewidth=0.8)
    bx.plot(distance, averaged, "o-", color=ACCENT, linewidth=1.8, markersize=5)
    bx.set_title("averaged over pairs at each separation", fontsize=9)
    bx.set_xlabel("separation $|i-j|$")
    bx.set_ylabel(r"$\langle S_i \cdot S_j \rangle$")
    bx.grid(True, color=MUTED, alpha=0.25, linewidth=0.6)

    _save(fig, "spin-correlations",
          {"matrix": spin.round(9).tolist(), "by_distance": averaged,
           "model": "hubbard L=8 U=4 (4,4)", "calculator": "alps_classic_dmrg",
           "max_bond_dim": 128})
    return {"nearest": averaged[0]}


def main() -> int:
    print("regenerating documentation figures")
    convergence()
    ci_wall()
    correlations()
    print("done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
