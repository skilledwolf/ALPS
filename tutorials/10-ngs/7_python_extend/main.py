 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # ALPS Project: Algorithms and Libraries for Physics Simulations                  #
 #                                                                                 #
 # ALPS Libraries                                                                  #
 #                                                                                 #
 # Copyright (C) 2010 - 2013 by Lukas Gamper <gamperl@gmail.com>                   #
 #                                                                                 #
 # ALPS Project: https://alps.comp-phys.org/                                       #
 # SPDX-License-Identifier: MIT                                                    #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import getopt
from pathlib import Path
import sys
import time
import pyalps.hdf5 as hdf5
import pyalps.ngs as ngs
import ising


if __name__ == '__main__':
    try:
        options, positional = getopt.getopt(sys.argv[1:], 'T:c')
        args = dict(options)
        limit = float(args.get('-T', 0))
        outfile = positional[0]
    except (IndexError, ValueError, getopt.GetoptError):
        sys.exit('usage: [-T timelimit] [-c] outputfile')

    simulation = ising.sim({'L': 100, 'THERMALIZATION': 100, 'SWEEPS': 1000, 'T': 2.})
    checkpoint = Path(outfile).with_suffix('.clone0.h5')
    clone = '/simulation/realizations/0/clones/0'
    if '-c' in args and checkpoint.exists():
        with hdf5.archive(checkpoint, 'r') as archive:
            archive.set_context(clone)
            simulation.load(archive)
    started = time.monotonic()
    simulation.run(lambda: limit > 0 and time.monotonic() > started + limit)

    def save(archive):
        archive.set_context(clone)
        simulation.save(archive)
    hdf5.save_checkpoint(str(checkpoint), save)
    results = ngs.collectResults(simulation)
    for name, result in results.items():
        print(f'{name}: {result.mean} +/- {result.error}')
    hdf5.save_checkpoint(outfile, lambda archive:
                         ngs.saveResults(results, simulation.parameters, archive, '/simulation/results'))
