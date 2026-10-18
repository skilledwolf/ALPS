# Copyright (c) 2026 Tobias Wolf
# SPDX-License-Identifier: MIT
#
# Adapted from skilledwolf/alpsx (private repository), branch
# agent/v3-execution-and-api-audit at d37e51227,
# python/qdk-chemistry-alps/examples/spin_quench.py. Changes: the classical
# trajectories are adaptive time-dependent DMRG from
# pyalps.classic_dmrg.evolve (second-order Trotter gates, kept states in
# place of an MPS bond cap) instead of alpsx's MPS evolution, with a
# state-fidelity check of its kept states; pyalps names (models.Hamiltonian,
# models.Spin, models.lattice.chain). The Q# part is unchanged.
"""One spin quench: exact propagation, ALPS tDMRG and executed Q# circuits."""

import time

import numpy as np
from demo_support import bit_reversal, finish, parser, pauli_matrix, site_matrix
from pyalps import classic_dmrg as dmrg
from pyalps import models
from pyalps.models import Hamiltonian, Spin
from qdk import qsharp

from qdk_chemistry_alps.interop import hamiltonian_to_qubit_operator

FIELDS = np.array([.07, -.13, .19, -.09, .16, -.21])
TIMES = np.arange(7) * .2


def model():
    spin = Spin()
    h = Hamiltonian(models.lattice.chain(6, bc="open"), spin)
    reference = np.zeros((64, 64), dtype=complex)
    for i in range(5):
        h.add(-2.8 * spin.Sz(i) * spin.Sz(i + 1))
        reference -= .7 * site_matrix(6, {i: "Z", i + 1: "Z"})
    for i, field in enumerate(FIELDS):
        h.add(-1.6 * spin.Sx(i) + 2 * field * spin.Sz(i))
        reference += -.8 * site_matrix(6, {i: "X"}) + field * site_matrix(6, {i: "Z"})
    return h, reference


def qsharp_source(operator):
    """Second-order symmetric product formula over the exported Pauli terms.

    Q# Exp(P, theta) implements exp(+i theta P); evolution therefore uses
    theta=-coefficient*dt/2. Pauli-container strings must be reversed when
    paired with a Q# array q[0], q[1], ... .
    """
    data = operator.to_dict()
    assert np.max(np.abs(data["coefficients"]["imag"])) < 1e-12
    gates = []
    for word, coefficient in zip(data["pauli_strings"], data["coefficients"]["real"], strict=True):
        paulis = ", ".join("Pauli" + p for p in word[::-1])
        gates.append(f"Exp([{paulis}], ({-float(coefficient):.17f}) * dt / 2.0, q);")
    return """operation QuenchEvolve(q : Qubit[], steps : Int, dt : Double) : Unit is Adj + Ctl {
    for _ in 1..steps {
        """ + "\n        ".join(gates + list(reversed(gates))) + """
    }
}
operation QuenchShot(steps : Int, dt : Double) : Result[] {
    use q = Qubit[6];
    QuenchEvolve(q, steps, dt);
    let results = [M(q[0]), M(q[1]), M(q[2]), M(q[3]), M(q[4]), M(q[5])];
    ResetAll(q);
    return results;
}
"""


def circuit_trajectory(source, dt, z_diagonal):
    qsharp.init()
    qsharp.eval(source)
    # StateDump indexes the Q# array in allocation order, q[0] MOST
    # significant, unlike QDK Chemistry's dense Pauli-register order.
    # Test this convention explicitly before any physics comparison.
    qsharp.eval("use q = Qubit[6]; X(q[0]);")
    probe = np.asarray(qsharp.dump_machine().as_dense_state())
    assert abs(probe[32] - 1.) < 1e-12
    qsharp.eval("X(q[0]);")
    values = []
    states = []
    steps = int(round(.2 / dt))
    for index, _ in enumerate(TIMES):
        if index:
            qsharp.eval(f"QuenchEvolve(q, {steps}, {dt:.17f});")
        state = np.asarray(qsharp.dump_machine().as_dense_state())
        np.testing.assert_allclose(np.linalg.norm(state), 1., atol=1e-11, rtol=0)
        values.append(abs(state)**2 @ z_diagonal)
        states.append(state)
    qsharp.eval("ResetAll(q);")
    return np.asarray(values), np.asarray(states)


def dmrg_trajectory(h, dt, kept):
    """Pauli Z of every site at TIMES, and the kept state of each, by tDMRG from |0...0>."""
    # Local state 0 of a spin 1/2 is up: Pauli Z = +1, Q# Zero.
    result = dmrg.evolve([0] * 6, h, time_step=dt, steps=round(TIMES[-1] / dt), states=kept,
                         measure_steps=round(.2 / dt), measurements={"Sz": "Sz"}, keep_mps=True)
    np.testing.assert_allclose(result.times, TIMES, atol=1e-12, rtol=0)
    np.testing.assert_allclose([state.norm() for state in result.mps], 1., atol=1e-10, rtol=0)
    return 2 * np.asarray(result.measurements["Sz"].values).real, np.array([s.to_dense() for s in result.mps])


def plot(plt, data):
    figure, axes = plt.subplots(1, 2, figsize=(11, 4), layout="constrained")
    image = axes[0].pcolormesh(TIMES, np.arange(6), np.asarray(data["exact"]).T, shading="nearest",
                               cmap="RdBu_r", vmin=-1, vmax=1)
    axes[0].set(xlabel="Time (ℏ/J₀)", ylabel="Site", xlim=(0, 1.2), title="Exact post-quench Pauli magnetization")
    figure.colorbar(image, ax=axes[0], label="⟨Zᵢ⟩ = 2⟨Sᶻᵢ⟩")
    for i, color in [(0, "#087e8b"), (2, "#d68c22"), (5, "#7956a8")]:
        axes[1].plot(TIMES, np.asarray(data["exact"])[:, i], color=color, label=f"Site {i}: exact")
        axes[1].plot(TIMES, np.asarray(data["dmrg_fine"])[:, i], "+", color=color)
        samples = np.asarray(data["shot_mean"])[:, i]
        low = np.asarray(data["shot_interval_low"])[:, i]
        high = np.asarray(data["shot_interval_high"])[:, i]
        axes[1].errorbar(data["sample_times"], samples, yerr=[samples - low, high - samples],
                         fmt="o", color=color, markersize=4, capsize=3)
    axes[1].set(xlabel="Time (ℏ/J₀)", ylabel="⟨Zᵢ⟩", ylim=(-1.1, 1.1),
                title="Lines: exact; +: tDMRG; points: Q# shots (95% CI)")
    axes[1].legend(fontsize=8)
    return figure


def main():
    arguments = parser(__doc__)
    arguments.add_argument("--shots", type=int, default=128, help="QDK shots per sampled time (default: 128)")
    args = arguments.parse_args()
    if args.shots < 1:
        arguments.error("--shots must be positive")
    started = time.perf_counter()
    h, exact_h = model()
    qdk = hamiltonian_to_qubit_operator(h)
    reverse = bit_reversal(6)
    np.testing.assert_allclose(pauli_matrix(qdk)[np.ix_(reverse, reverse)], exact_h, atol=1e-12, rtol=0)
    for i in range(6):
        observable = Hamiltonian(h.lattice, h.site_for(i)).add(2 * h.site_for(i).Sz(i))
        mapped = pauli_matrix(hamiltonian_to_qubit_operator(observable))
        np.testing.assert_allclose(mapped[np.ix_(reverse, reverse)], site_matrix(6, {i: "Z"}), atol=1e-12, rtol=0)
    eigenvalues, vectors = np.linalg.eigh(exact_h)
    states = (np.exp(-1j * TIMES[:, None] * eigenvalues) * vectors[0].conj()) @ vectors.T
    z_diagonal = np.column_stack([np.diag(site_matrix(6, {i: "Z"})).real for i in range(6)])
    exact = abs(states)**2 @ z_diagonal
    source = qsharp_source(qdk)
    circuit_coarse, _ = circuit_trajectory(source, .1, z_diagonal)
    circuit_fine, circuit_states = circuit_trajectory(source, .05, z_diagonal)
    classical_coarse, _ = dmrg_trajectory(h, .04, 8)
    classical_fine, classical_states = dmrg_trajectory(h, .02, 8)
    classical_few_states, _ = dmrg_trajectory(h, .02, 2)
    errors = {"qdk_dt_0.1": np.max(abs(circuit_coarse - exact)),
              "qdk_dt_0.05": np.max(abs(circuit_fine - exact)),
              "dmrg_dt_0.04_states_8": np.max(abs(classical_coarse - exact)),
              "dmrg_dt_0.02_states_8": np.max(abs(classical_fine - exact)),
              "dmrg_dt_0.02_states_2": np.max(abs(classical_few_states - exact))}
    assert errors["qdk_dt_0.05"] < errors["qdk_dt_0.1"]
    assert errors["dmrg_dt_0.02_states_8"] < errors["dmrg_dt_0.04_states_8"]
    assert errors["dmrg_dt_0.02_states_8"] < errors["dmrg_dt_0.02_states_2"]
    np.testing.assert_allclose(circuit_fine, exact, atol=.004, rtol=0)
    np.testing.assert_allclose(classical_fine, exact, atol=.001, rtol=0)
    # A phase-sensitive state check catches the wrong sign of Exp, even when
    # real-Hamiltonian Z trajectories alone cannot distinguish t from -t.
    overlaps = np.abs(np.sum(states.conj() * circuit_states, axis=1))**2
    assert np.min(overlaps) > .9999
    # The kept tDMRG states are site 0 most significant, as the reference.
    classical_overlaps = np.abs(np.sum(states.conj() * classical_states, axis=1))**2
    assert np.min(classical_overlaps) > 1 - 1e-6

    sample_indices = [2, 4, 6]
    samples = []
    for index in sample_indices:
        results = qsharp.run(f"QuenchShot({round(TIMES[index] / .05)}, 0.05)",
                             shots=args.shots, seed=421 + index)
        if any(not isinstance(row, list) or len(row) != 6 for row in results):
            raise RuntimeError(f"QDK shot failed: {results!r}")
        samples.append([[1 if bit == qsharp.Result.Zero else -1 for bit in row] for row in results])
    samples = np.asarray(samples)
    means = samples.mean(axis=1)
    # Family-wise Hoeffding bound at failure probability 0.001 for 18 means.
    bound = np.sqrt(2 * np.log(2 * means.size / .001) / args.shots)
    assert np.max(abs(means - circuit_fine[sample_indices])) <= bound
    # Wilson intervals stay meaningful at all-equal outcomes.
    p, n, z = (means + 1) / 2, args.shots, 1.959963984540054
    center = (p + z*z / (2*n)) / (1 + z*z/n)
    radius = z * np.sqrt(p*(1-p)/n + z*z/(4*n*n)) / (1 + z*z/n)
    print("Six-site Ising quench, actual Q# Exp circuits and Z measurements")
    print("Maximum magnetization errors against independent dense propagation:",
          ", ".join(f"{name} {value:.2g}" for name, value in errors.items()))
    print(f"Minimum QDK state fidelity: {overlaps.min():.9f}; tDMRG: {classical_overlaps.min():.9f}")
    print(f"{args.shots} shots × 3 times, 95% Wilson intervals")
    finish(args, {"model": {"sites": 6, "ZZ": -.7, "X": -.8, "Z": FIELDS, "boundary": "open"},
                  "initial_state": "all ALPS local indices 0 = all Q# qubits Zero = Pauli Z +1",
                  "units": "hbar=1; Hamiltonian coefficients in J0; time in hbar/J0",
                  "times": TIMES, "exact": exact, "dmrg_fine": classical_fine,
                  "dmrg_coarse": classical_coarse, "dmrg_few_states": classical_few_states,
                  "qdk_fine": circuit_fine, "qdk_coarse": circuit_coarse, "errors": errors,
                  "state_fidelity": overlaps, "dmrg_state_fidelity": classical_overlaps, "qsharp_source": source,
                  "controls": {"dmrg_dt": [.04, .02], "dmrg_states": [2, 8], "dmrg_trotter_order": 2,
                               "qdk_dt": [.1, .05], "qdk_order": 2, "shots": args.shots, "seeds": [423, 425, 427]},
                  "sample_times": TIMES[sample_indices], "shot_mean": means,
                  "shot_plus_counts": (samples == 1).sum(axis=1),
                  "shot_interval_low": 2*(center-radius)-1, "shot_interval_high": 2*(center+radius)-1,
                  "familywise_shot_error_bound": bound}, started, plot)


if __name__ == "__main__":
    main()
