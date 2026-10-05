# Spatial heat-bath Ising

`ising_multiple` uses native ALEA and the same TOML schema, statistics,
temperature scans and restart rules as [ising_single](../single/README.md).
Build it from the parent examples project against an installed ALPS SDK:

```sh
build/parapack/multiple/ising_multiple --validate parapack/multiple/run.toml
mpiexec -n 3 build/parapack/multiple/ising_multiple parapack/multiple/run.toml
```

With an MPI-enabled SDK, all ranks cooperate on each physical chain. A periodic
nearest-neighbor chain is divided into contiguous blocks with at least two
sites per rank. Odd lengths and uneven partitions are supported. Ghost spins
are exchanged after each color update; extensive energy and magnetization are
reduced over the blocks. `execution.chains` selects independent chains, each
using the whole communicator in turn. Spatial ranks do not multiply the
number of samples. Other graphs are supported when running on one rank or
using a serial SDK.

Optional OpenMP updates within each block follow the same fixed proposal order.
Enable `ALPS_ENABLE_OPENMP` for the examples and use `OMP_NUM_THREADS`; MPI calls
stay on the main thread. Thread counts also preserve results and checkpoints.

Proposals follow a fixed global order. Native checkpoints contain canonical
global spins, ordered topology, the RNG and all partially accumulated statistics.
Restart with a different valid rank count preserves the trajectory and results,
including scans and either supported RNG. All ranks agree on time/signal stopping
before each sweep, gather spins for checkpoints, and publish through rank zero.

The command retains the ring Hamiltonian `-J*sum_b s_i*s_j`, with spins ±1 and
either sign of `J`. The historical scheduler output goldens used incorrect
acceptance weights and energy signs; they have been removed. The corrected
legacy worker in `ising.h` remains only for the unported nested replica-exchange
adapter. Released Parapack physical checkpoints still need an offline converter.
