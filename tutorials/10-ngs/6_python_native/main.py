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

import pyalps.hdf5 as hdf5
import pyalps
import getopt
from pathlib import Path
import sys
import time

import ising

if __name__ == '__main__':

    try:
        optlist, positional = getopt.getopt(sys.argv[1:], 'T:c')
        args = dict(optlist)
        limit = float(args.get('-T', 0))
        resume = '-c' in args
        outfile = positional[0]
    except (IndexError, ValueError, getopt.GetoptError):
        sys.exit('usage: [-T timelimit] [-c] outputfile')

    sim = ising.sim({
        'L': 100,
        'THERMALIZATION': 100,
        'SWEEPS': 1000,
        'T': 2
    })

    checkpoint = Path(outfile).with_suffix('.clone0.h5')
    if resume and checkpoint.exists():
        with hdf5.archive(checkpoint, 'r') as ar:
            sim.load(ar)

    if limit == 0:
        sim.run(lambda: False)
    else:
        start = time.monotonic()
        sim.run(lambda: time.monotonic() > start + limit)

    hdf5.save_checkpoint(str(checkpoint), sim.save)

    results = sim.collectResults()
    for key, value in results.items():
        print(f'{key}: {value.mean} +/- {value.error}')

    def save_results(ar):
        ar['/parameters'] = sim.parameters
        for name, value in results.items():
            ar['/simulation/results/' + pyalps.hdf5_name_encode(name)] = value

    hdf5.save_checkpoint(outfile, save_results)
