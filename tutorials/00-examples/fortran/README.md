# Fortran Monte Carlo examples

These examples use the installed `ALPS::fortran` target, typed TOML parameters,
`alps::mcbase`, native ALEA batches and autocorrelation analysis. A Fortran compiler
supporting `ISO_C_BINDING` is required for the examples; the SDK bridge itself
only needs C++. No OpenMP or scheduler is required. The optional tests require
Python 3.11 or later with NumPy and h5py.

Build this directory against an installed ALPS package. Run `hello_fortran
hello/run.toml` or `ising_fortran ising/run.toml`; `--schema` and `--validate` use
the same CLI as the native MC applications. The corresponding ALPSize lessons
build these same sources under their original command names, `hello` and `ising`.
Output paths are relative to the run file. Set `execution.chains` for independent
chains and `execution.max_sweeps` for a bounded invocation. Resume with
`input.checkpoint` and distinct output paths. Increasing `parameters.SWEEPS`
extends production while retaining the saved RNG, samples and physical state.

The Ising example implements sequential heat-bath updates on a periodic square
lattice, with Hamiltonian `-sum_<ij> s_i s_j`. Energy and magnetization are per
site. `THERMALIZATION` updates precede `SWEEPS` production updates. This replaces
the old `INT`/`MCS` parameters and overflow-prone example RNG; trajectories change.
The energy sign is corrected to agree with the Hamiltonian.

The hello example also demonstrates vector observations, character arrays,
32- and 64-bit integers, and single- and double-precision checkpoint data.

## Writing an application

Use `iso_c_binding`, then include `alps/fortran/alps_fortran.h` after `implicit
none`. Callbacks use `bind(C)` and `type(c_ptr), value :: caller`. All names passed
to the bridge are NUL-terminated (`'Energy'//c_null_char`). Data pointers use
`c_loc` and the corresponding interoperable type: `c_int`, `c_int64_t`, `c_float`,
`c_double` or `c_char`. `ALPS_INT64` denotes a fixed 64-bit integer on all platforms.
Counts and character widths use `c_size_t`; numeric widths are zero. Character
buffers and arrays use an explicit width, with blank padding on reads.

Allocate one state object in `alps_init`, attach it using `alps_set_context`, and
retrieve it with `alps_get_context`/`c_f_pointer`. `alps_finalize` must tolerate a
partly initialized state. Do not keep mutable state in module globals: several
chains and staged checkpoint loads can coexist in one process. Draw randomness
through `alps_random(caller)`.

Register scalar or vector observables in `alps_init_observables`. Accumulate them
in `alps_run`; the bridge discards observations when `alps_is_thermalized` was
false **before** that update. Different observables may use different sampling
intervals. Pooled output preserves each chain's weighted batches; chronological
autocorrelation data remain separate for each chain.

Save physical state with ordered `alps_dump` calls and load it with matching
`alps_restore` calls. HDF5 stores typed arrays, a version and field descriptors;
there is no private binary encoding. The framework saves the RNG and native
accumulators. Validate application invariants after loading and report failures
using `alps_fail`. Check `alps_failed(caller)` and return before using failed
reads: C++ exceptions are retained and rethrown after the callback returns.
Loads are staged in a separate instance, so a rejected checkpoint leaves the
live simulation intact. The caller handle itself must never be stored in state.

The old F77 callback ABI, scheduler input files and XDR physical dumps are not
accepted by this bridge. Existing recorded HDF5 results can use the offline
archive converter. Old application-defined, untyped XDR state needs an
application-specific offline decoder; generic native restart from it is not
implemented. Keep those source checkpoints until that conversion is available.
