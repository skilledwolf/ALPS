# Replica exchange with spatially decomposed walkers

`exchange` runs a replica-exchange ladder of heat-bath Ising walkers with
uniform coupling `J`. It accepts the parameters of the classical `ising`
example. Each walker can also be spread over several MPI processes, which
the classical examples cannot do. With the installation's `bin` directory on
`PATH`:

```sh
exchange --validate run.toml
exchange run.toml
mpiexec -n 4 exchange spatial.toml
```

`run.toml` follows the former `params-ising` job: a 16-temperature ladder on
an 8×8 square lattice with temperature optimization. `spatial.toml`
decomposes each 64-site ring walker over two processes, so four processes
run two walkers at a time.

## Execution layouts

| `execution` settings | Processes | Layout |
| --- | --- | --- |
| `parallel = "chains"` (default) | any | Independent ladders are distributed over processes |
| `parallel = "replicas"` | any | The walkers of one ladder are distributed over processes |
| `parallel = "replicas"`, `processes_per_walker = p` | a multiple of `p` | Consecutive groups of `p` processes form teams; walker `w` runs on team `w % teams`, spatially decomposed over its `p` processes |

Spatial decomposition needs a periodic nearest-neighbor ring (`chain
lattice`) with at least two sites per process. Other lattices run with
`processes_per_walker = 1`. Each walker draws all of its proposals from one
random stream, so results and checkpoints are identical for every layout;
a checkpoint can resume with a different layout and process count.

## Former algorithms

| Former `ALGORITHM` | Native equivalent |
| --- | --- |
| `"ising; exchange"` | The same name, any layout above |
| `"multiple parallel ising; exchange"` with `PROCESS_PER_WORKER = p` | The same name with `execution.processes_per_walker = p` |
| `"loop; exchange"` | The `loop` application, which supports replica exchange and `parallel = "replicas"` |

Observables and moments are those of `ising_single`: energy and magnetization
are extensive. Each temperature's results are under
`/simulation/replicas/<id>/results`, with the exchange diagnostics of the
[classical examples](../../README.md#native-classical-examples).
Released Parapack checkpoints are not converted; finish those runs with
ALPS 3.0.
