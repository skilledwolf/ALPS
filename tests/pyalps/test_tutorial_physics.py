"""Check actual tutorial updates against independently calculated physics."""
import importlib.util
import itertools
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest

TUTORIALS = Path(__file__).resolve().parents[2] / "tutorials"


def load_tutorial(relative):
    path = TUTORIALS / relative
    spec = importlib.util.spec_from_file_location("tutorial", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize("relative", [
    "ngs/6_python_native/ising.py", "ngs/7_python_extend/ising.py",
    "code-08-mcmain-python/ising.py", "code-09-mcmain-python-hybrid/ising.py",
])
def test_ising_single_flip_obeys_boltzmann_weights(relative):
    tutorial = load_tutorial(relative)

    class OneProposalFinished(Exception):
        pass

    # Use the actual update method with a controllable RNG. Abort at the next
    # proposal, so every state/site and both sides of its acceptance threshold
    # can be checked without Monte Carlo noise or native RNG assumptions.
    class State:
        length = 4
        beta = 0.7

    def energy(spins):
        return -np.dot(spins, np.roll(spins, 1))

    for spins in itertools.product([-1., 1.], repeat=4):
        before = np.array(spins)
        for site in range(4):
            after = before.copy()
            after[site] *= -1
            delta = energy(after) - energy(before)
            acceptance = min(1., np.exp(-State.beta * delta))
            for draw in [acceptance / 2, (1 + acceptance) / 2]:
                state = State()
                state.spins = before.copy()
                draws = iter([(site + .25) / 4] + ([draw] if acceptance < 1 else []))
                def random():
                    try:
                        return next(draws)
                    except StopIteration:
                        raise OneProposalFinished
                state.random = random
                with pytest.raises(OneProposalFinished):
                    tutorial.sim.update(state)
                np.testing.assert_array_equal(
                    state.spins, after if acceptance == 1 or draw < acceptance else before,
                    err_msg=f"{relative}, spins={spins}, site={site}, draw={draw}")


def test_density_profile_runs_and_plots_in_scratch_directory(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    tutorial = load_tutorial("optical-lattice-02-density-profile/density_profile.py")
    assert not list(tmp_path.iterdir())  # Importing must not launch a simulation.
    parms = tutorial.parameters()
    assert parms["L"] == 9 and parms["MEASURE_LOCAL[Local Density]"] == "n"
    trap = parms["mu"]
    assert eval(trap, {"__builtins__": {}}, dict(x=4, y=4, z=4)) == pytest.approx(4.05)
    assert eval(trap, {"__builtins__": {}}, dict(x=0, y=0, z=0)) == pytest.approx(4.05 - .65 * 48)
    # At Nmax=1, mu=0, particle-hole symmetry fixes <n_i>=1/2 at every
    # site of the homogeneous hopping model, independent of temperature.
    # EPSILON supplies the diagonal offset required for purely hopping SSE.
    # A fixed seed and generous absolute tolerance make this an integration
    # regression, not a convergence or uncertainty-estimator test.
    density = tutorial.run(length=2, curvature=0., Nmax=1, mu=0.,
                           THERMALIZATION=500, SWEEPS=10000, SEED=42, EPSILON=.5)
    assert density.shape == (2, 2, 2)
    np.testing.assert_allclose(density, .5, rtol=0., atol=.08)
    fig = tutorial.plot(density)
    np.testing.assert_allclose(fig.axes[0].lines[0].get_ydata(), density[:, 1, 1])
    np.testing.assert_allclose(fig.axes[1].images[0].get_array(), density[:, :, 1])
    fig.savefig(tmp_path / "density.png")
    assert (tmp_path / "density.png").stat().st_size > 0
    plt.close(fig)
