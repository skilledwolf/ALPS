The three examples use `alps::mcbase` and native ALEA batch measurements:

- [A periodic chain with three-component spins](1d_lattice/).
- [Three-component spins on an ALPS lattice](nd_lattice/).
- [An arbitrary number of spin components on an ALPS lattice](o_n_model/).

Build each directory with CMake against an installed ALPS SDK. Run the executable with its TOML parameter file; the `test/param.toml` files provide lattice examples. Results and resumable batches use canonical HDF5 serialization.

The ferromagnetic energy is `-sum(dot(spin_i, spin_j))`. A proposed spin therefore has acceptance probability `min(1, exp(beta * dot(new_spin - old_spin, neighbors)))`.

For volume `V`, susceptibility is the ensemble fluctuation `beta * V * (<M²> - dot(<M>, <M>))`. Compute it from the recorded magnetization and its second moment rather than an instantaneous subtraction. The arbitrary-component example also records distances and correlations from site zero.
