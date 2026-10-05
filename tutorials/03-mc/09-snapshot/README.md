# Classical spin snapshots

Install pyalps and the native `simplemc` executable, then run `python run.py b`
for Ising snapshots, `python run.py c` for XY or `python run.py d` for
Heisenberg. Each writes explicit TOML tasks and native HDF5 results, plus VTK
snapshots that can be opened directly in ParaView. The included `plot9b.pvsm`,
`plot9c.pvsm` and `plot9d.pvsm` state files use those VTK filenames.

`python run.py a` runs the finite-size temperature sweep; `python plot9a.py`
plots its specific heat, magnetization square and Binder ratio. `run.sh` runs all
four groups. `clean.sh` removes their generated task, result and snapshot files.

The [simplemc guide](../../../src/apps/mc/simple/README.md) explains native
statistics, independent chains and checkpoint restart. VTK is emitted directly;
`snap2vtk` is only needed to convert old released `.snap` files.
