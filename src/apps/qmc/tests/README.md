# QMC settings contracts

These tests characterize each solver's existing input and checkpoint rules.
The settings extraction preserves those differences for drivers that need to
resolve parameters without constructing a lattice or a scheduler worker.
Looper's `temperature` and `mc_steps` already provide that separation.

| Solver | Temperature keys, in precedence order | Sweep count | Thermalization | Changing settings on resume |
| --- | --- | --- | --- | --- |
| worm | `Beta`, `beta`, `BETA`, `T`, `TEMPERATURE`, `temperature`; expressions evaluated | `SWEEPS`, required plain number at construction; negative or non-finite rejected | Required plain integer; negative rejected | Worker rejects changed values except `SEED`; `SWEEPS` cannot be extended |
| dirloop_sse | Same temperature rule as worm | Plain unsigned integer: `SWEEPS`, `MCS`, `Steps` | `THERMALIZATION`, `thermalization`, then integer `sweeps / 10` | Positive sweep count changes accepted; checkpoint restores its saved thermalization |
| loop | `T`, `TEMPERATURE`, `BETA`, `INVERSE_TEMPERATURE`; expressions evaluated; lowercase `beta` ignored | Integer range, default `[65536:]` | Expression if present, otherwise minimum sweeps divided by eight | `mc_steps` checkpoint restores its schedule; driver can call `set_sweeps` to extend it |
| qwl | Ignored; beta is fixed at 1 | Wang-Landau schedule, then multicanonical `SWEEPS` | No thermalization setting | Application accepts only a changed `SWEEPS` (worker also permits `SEED`) |

For comparison, spinmc's temperature rule is lowercase `beta` before `T`, with
plain numeric conversion. This is intentionally different from all four
contracts above. Looper supports annealing segments in its temperature helper.
The `mc_steps` helper wraps negative thermalization to unsigned; the loop
workers then reject it through their annealing-length check because its
`thermalization()` getter returns a signed integer. Negative SSE counts also
wrap. These are existing rules, not validation changes introduced by the
extraction.

Run the tests from a build configured with `ALPS_BUILD_TESTING=ON` and
`ALPS_BUILD_APPLICATIONS=ON`:

```sh
cmake --build build --target qmc_worm_characterization qmc_sse_characterization \
  qmc_qwl_characterization qmc_loop_characterization \
  qmc_worm_settings qmc_sse_settings qmc_qwl_settings
ctest --test-dir build -R '^qmc\.' --output-on-failure
```

The assertions cover settings, progress, and checkpoint restoration, without
committed reference measurements that would couple these tests to independent
physics corrections. Optional fixed-step dumps also record RNG, measurements,
and solver state to compare two builds on the same platform:

```sh
ALPS_QMC_DUMP_DIR=/tmp/qmc-before ctest --test-dir build -R '^qmc\.' --output-on-failure
# Save the unchanged binaries, then rebuild after the extraction.
ALPS_QMC_DUMP_DIR=/tmp/qmc-after ctest --test-dir build -R '^qmc\.' --output-on-failure
python3 src/apps/qmc/tests/compare_dumps.py /tmp/qmc-before /tmp/qmc-after
```

The optional comparison requires `h5py` and `numpy`. Run the unchanged build
twice first to confirm repeatability. Numeric HDF5 datasets are compared by
their value bytes; strings by content; group and dataset names and all
attributes are also compared; SSE and loop solver XDR payloads by bytes.
Worker log timestamps are excluded by saving solver payloads directly.

The scheduler accepts a changed `SEED` parameter without reseeding a resumed
worker; the application's `change_parameter` hook is not called for it.

Two preexisting uninitialized fields prevent comparing every raw checkpoint
byte. Worm's `last_id_` counter is also copied into kink IDs. Its canonical
HDF5 snapshot retains kink times, states, site counts, and IDs relative to
`last_id_`, preserving pairing and relative ID history. QWL initializes only
the opcode of identity operators; the unused bond of each identity operator is
normalized to zero in its canonical HDF5 snapshot. Active operator bonds are
retained. All other decoded state is compared without normalization.

QWL can be constructed with a positive `EXPANSION_ORDER_MINIMUM`, but stepping
that configuration crashes in the unchanged solver. The tests characterize
construction only for that case. Fixing it is separate from preserving settings.
