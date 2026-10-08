# Command-line tools

These directories own executable consumers of ALPS modules. They do not add
new independently linkable libraries. The root CMake file dispatches to each
group; executable names and installation components remain stable.

| Group | Built or installed programs | Preserved inactive sources |
| --- | --- | --- |
| `lattice/` | `lattice2xml`, `printgraph` | `pltgraph.py` |
| `snapshot/` | `snap2vtk` | None |
| `diagnostics/` | `pconfig` | None |
| `hdf5/` | `alps-hdf5-convert` (Python, requires h5py and NumPy) | None |
| `xml/` | `alps-xml` on Unix | `txt2archive.C` and all historical shell wrappers |
| `alea/` | None | C++ and Python mean/variance analysis programs |
| `launchers/` | None | Historical `alpspython` shell and Windows templates |
| `installer/` | None | Historical macOS postflight template |

All four C++ commands and `alps-hdf5-convert` install in the `tools` component;
`alps-xml` installs in the `xml` component. `pconfig` links `ALPS::utilities`; the other C++ commands
continue to use `ALPS::alps`. Source ownership is registered for all C++ groups,
including inactive files, without enabling additional programs or dependencies.

`printgraph run.toml` and `lattice2xml run.toml` read the lattice from the
`[parameters]` of a TOML run file, with an optional `input.lattice_library`, as the
lattice-model applications do; `lattice2xml LATTICE [KEY=VALUE ...]` remains as a
short form. The result archive group is separate from the HDF5
implementation. Inactive sources are preserved for a separate functionality
assessment; their presence does not claim current build or runtime support.
Central CLI integration tests and their fixtures remain under `tests/cli/`.

The [HDF5 converter](hdf5/README.md) is independent of the ALPS runtime. It moves
legacy complex/Boolean encodings to ordinary HDF5 datatypes in a separate file,
isolating migration from the future HighFive implementation.
