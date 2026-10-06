# Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
"""Run the temperature scan or independent depleted-lattice realizations."""
import argparse

from pyalps.run_io import execute, write_run_files


def runs(mode):
    base = dict(ALGORITHM='loop', MODEL='spin', local_S=.5, J=1.)
    if mode == 'scan':
        parameters = [dict(base, LATTICE='chain lattice', L=8, T=t,
                           SWEEPS=65536, THERMALIZATION=6553)
                      for t in (.1, .2, .3, .4, .5)]
    else:
        parameters = [dict(base, LATTICE='depleted square lattice', L=4, T=1.,
                           DEPLETION=.2, DEPLETION_SEED=23123 + i,
                           SWEEPS=128, THERMALIZATION=12)
                      for i in range(100)]
    return [dict(parameters=p,
                 execution=dict(seed=42 + i, disorder_seed=23123 + i, chains=1),
                 output=dict(results=f'{mode}{i}.h5', checkpoint=f'{mode}{i}.checkpoint.h5'))
            for i, p in enumerate(parameters)]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=['scan', 'disorder'])
    mode = parser.parse_args().mode
    execute('loop', write_run_files(mode, runs(mode)))
