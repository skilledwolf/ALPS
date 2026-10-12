# ALPS Project: https://alps.comp-phys.org/
# SPDX-License-Identifier: MIT
"""Density profile of trapped bosons in a cubic optical lattice.

Simulates the Bose-Hubbard model with a parabolic trap,

    H = -t sum_<ij> b_i^+ b_j + U/2 sum_i n_i (n_i - 1) - sum_i (mu - V_T(r_i)) n_i,

with the directed-loop SSE quantum Monte Carlo code and plots the local density
<n_i>. U/t = 8.11 is the value optical-lattice-01-bandstructure gives for a
V0 = 8 E_r lattice of 87Rb. This replaces the dwa-02-density-profile tutorial,
which was removed together with the DWA application.

The expected output is shown in density_profile.png. The lattice is kept small
so the tutorial runs in a few minutes, which means it still sits inside the
cloud: the density falls off toward the edges but does not quite vanish
there. A bigger lattice (larger L) lets the density drop closer to zero at
the edges.
"""

import numpy as np
import matplotlib.pyplot as plt
import pyalps

def parameters(length=9, curvature=0.65, **overrides):
    """Construct the trapped Bose-Hubbard inputs; overrides allow small runs."""
    center = (length - 1) / 2.
    result = {
        'LATTICE': 'inhomogeneous simple cubic lattice',
        'L': length,
        'MODEL': 'boson Hubbard',
        'Nmax': 4,
        't': 1.,
        'U': 8.11,
        'mu': '4.05 - %g*((x-%g)*(x-%g) + (y-%g)*(y-%g) + (z-%g)*(z-%g))'
              % ((curvature,) + (center,) * 6),
        'T': 1.,
        'THERMALIZATION': 1000,
        'SWEEPS': 5000,
        'MEASURE_LOCAL[Local Density]': 'n',
    }
    result.update(overrides)
    return result


def run(length=9, curvature=0.65, prefix='parm_trap', **overrides):
    """Run the simulation in the working directory and return its density cube."""
    parms = parameters(length, curvature, **overrides)
    input_file = pyalps.writeInputFiles(prefix, [parms])
    status, _ = pyalps.runApplication('dirloop_sse', input_file, Tmin=1)
    if status != 0:
        raise RuntimeError('dirloop_sse failed with status %s' % status)
    data = pyalps.loadMeasurements(pyalps.getResultFiles(prefix=prefix), 'Local Density')[0][0]
    return np.asarray(data.y.mean).reshape(length, length, length)


def plot(density):
    """Plot the central cut and layer without opening an interactive window."""
    length = density.shape[0]
    r = np.arange(length) - (length - 1) / 2.
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
    ax1.plot(r, density[:, length // 2, length // 2], 'o-')
    ax1.set_xlabel('$x$ (lattice sites from trap center)')
    ax1.set_ylabel(r'$\langle n_i \rangle$')
    ax1.set_title('Cut through the trap center')
    im = ax2.imshow(density[:, :, length // 2], extent=[r[0], r[-1], r[0], r[-1]], origin='lower')
    fig.colorbar(im, ax=ax2, label=r'$\langle n_i \rangle$')
    ax2.set_title('Center layer $z = 0$')
    fig.tight_layout()
    return fig


def main():
    density = run()
    center = density.shape[0] // 2
    print('Total number of bosons: %.2f' % density.sum())
    print('Density at the trap center: %.3f' % density[center, center, center])
    plot(density)
    plt.show()


if __name__ == '__main__':
    main()
