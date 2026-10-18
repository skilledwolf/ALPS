<!--
Copyright (c) 2026 Tobias Wolf
SPDX-License-Identifier: MIT
-->

# Figures

The figures are committed, so reading the documentation needs neither matplotlib nor the solver stack. Each figure that shows computed values sits next to the JSON file it was plotted from.

| Figure | Page | Made by |
| --- | --- | --- |
| `dmrg-convergence.svg`, `.json` | [Converging DMRG](../how-to/converging-dmrg.md) | `tools/generate_figures.py` |
| `ci-wall.svg`, `.json` | [Large active spaces](../how-to/large-active-spaces.md) | `tools/generate_figures.py` |
| `spin-correlations.svg`, `.json` | [Observables](../how-to/observables.md) | `tools/generate_figures.py` |
| `mode-ordering.svg` | [Conventions](../conventions.md) | Drawn by hand; copied from alpsx with a white background added (provenance in the file) |
| `demos/*.png`, `demos/*.json` | [Examples](../examples/index.md) | The example scripts, with `--output` |

## Regenerate

After [installing the solver stack](../install.md) and matplotlib, run from the **ALPS repository root**:

```bash
python python/qdk-chemistry-alps/tools/generate_figures.py
```

This rewrites the three data figures in place. The SVGs carry no date, so an unchanged result gives an unchanged file.

The gallery figures come from the examples:

```bash
python python/qdk-chemistry-alps/examples/hubbard_reference.py --output DIRECTORY
```

and likewise `custom_spin_export.py` and `fermion_pairing.py`. Each run writes `<script>.png` and `<script>.json`; the JSON records the model, the numerical controls, the dependency versions, the git revision and whether the tree was dirty. Use a directory outside the repository and copy the files into `demos/`: a second run into `demos/` would see the first run's files as untracked and record a dirty tree.
