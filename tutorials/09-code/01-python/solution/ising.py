# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2009-2010 by Matthias Troyer <troyer@phys.ethz.ch> 
#                            Jan Gukelberger
#                            Brigitte Surer
# 
# ALPS Project: https://alps.comp-phys.org/
# SPDX-License-Identifier: MIT
# 
# ****************************************************************************

import argparse
import math
import operator
from pathlib import Path

import numpy as np
import pyalps
from pyalps import alea, hdf5, run_config
from pyalps.ngs import params, random01

SCHEMA = Path(__file__).with_name('schema.toml').read_text()
NAMES = ('E', 'm', '|m|', 'm^2', 'm^4')


class Simulation:
    def __init__(self, beta, L, seed=42, bins=64):
        run = run_config.resolve(SCHEMA, parameters={'BETA': beta, 'L': L},
                                 execution={'seed': seed, 'bins': bins},
                                 output={'results': 'unused.h5'})
        self.parameters = params(dict(run.parameters))
        self.L, self.beta = self.parameters['L'], self.parameters['BETA']
        self.rng = random01(run.execution['seed'])
        # Only uphill moves need an exponential: avoid overflow at large beta.
        self.exp_table = {e: math.exp(2*self.beta*e) for e in (-4, -2, 0)}
        self.spins = [[2*self.randint(2)-1 for _ in range(self.L)] for _ in range(self.L)]
        # Joint samples retain correlations needed by nonlinear observables.
        self.samples = alea.BatchAccumulator(len(NAMES), num_batches=run.execution['bins'])
        self.diagnostics = {name: alea.AutocorrelationAccumulator() for name in NAMES}
        self.execution = run.execution
        self.started = False

    def run(self, ntherm, n):
        ntherm, n = operator.index(ntherm), operator.index(n)
        if self.started or ntherm < 0 or n < 2:
            raise ValueError('Use a fresh simulation, nonnegative warmup and at least two production sweeps')
        self.started = True
        self.parameters['THERMALIZATION'], self.parameters['SWEEPS'] = ntherm, n
        for _ in range(ntherm):
            self.step()
        for _ in range(n):
            self.step()
            sample = self.observables()
            self.samples << sample
            for name, value in zip(NAMES, sample):
                self.diagnostics[name] << value
        result = self.samples.result()
        for i, name in enumerate(NAMES[:3]):
            diagnostics = self.diagnostics[name].result()
            tau = diagnostics.tau[0] if diagnostics.tau_available else 'unavailable (insufficient bins)'
            print(f'{name}: {result.mean[i]} +/- {result.error[i]}, tau = {tau}')

    def step(self):
        for _ in range(self.L*self.L):
            i, j = self.randint(self.L), self.randint(self.L)
            neighbors = (self.spins[(i-1)%self.L][j] + self.spins[(i+1)%self.L][j]
                         + self.spins[i][(j-1)%self.L] + self.spins[i][(j+1)%self.L])
            e = -self.spins[i][j]*neighbors
            if e > 0 or self.rng() < self.exp_table[e]:
                self.spins[i][j] = -self.spins[i][j]

    def observables(self):
        energy = magnetization = 0.
        for i in range(self.L):
            for j in range(self.L):
                energy -= self.spins[i][j]*(self.spins[(i+1)%self.L][j] + self.spins[i][(j+1)%self.L])
                magnetization += self.spins[i][j]
        e, m = energy/self.L**2, magnetization/self.L**2
        return np.array([e, m, abs(m), m*m, m**4])

    def randint(self, maximum):
        return int(maximum*self.rng())

    def results(self):
        joint = self.samples.result()
        results = {name: joint.transform(lambda x, i=i: x[i:i+1]) for i, name in enumerate(NAMES)}
        def ratio(x):
            if x[3] <= 0:
                raise ValueError('Binder ratio requires positive <m^2> in every jackknife sample')
            return np.array([x[4]/x[3]**2])
        unavailable = {}
        try:
            results['Binder Ratio'] = joint.transform(ratio)
        except ValueError as error:
            unavailable['Binder Ratio'] = str(error)
        return results, unavailable

    def save(self, filename):
        if not self.started or self.samples.count != self.parameters['SWEEPS']:
            raise ValueError('Complete the run before saving its analysis results')
        results, unavailable = self.results()
        run = run_config.resolve(SCHEMA, parameters=dict(self.parameters),
                                 execution=dict(self.execution), output={'results': str(filename)})
        def write(ar):
            ar['/parameters'] = self.parameters
            ar['/run_config'] = run
            for name, result in results.items():
                result.save(ar, '/simulation/results/' + pyalps.hdf5_name_encode(name))
            for name, reason in unavailable.items():
                ar['/simulation/unavailable/' + pyalps.hdf5_name_encode(name)] = reason
            self.samples.result().save(ar, '/simulation/joint')
            for name, accumulator in self.diagnostics.items():
                accumulator.result().save(ar, '/simulation/realizations/0/clones/0/autocorrelation/'
                                          + pyalps.hdf5_name_encode(name))
        hdf5.save_checkpoint(str(filename), write)


def main(simulation=Simulation, sizes=(4,), plot=False):
    parser = argparse.ArgumentParser(description='Square-lattice Ising Monte Carlo with native ALEA')
    parser.add_argument('--schema', action='store_true')
    parser.add_argument('--validate', action='store_true')
    parser.add_argument('runs', nargs='*', help='TOML run files; no arguments runs the tutorial scan')
    args = parser.parse_args()
    if args.schema:
        print(SCHEMA)
        return
    if args.validate and not args.runs:
        parser.error('--validate requires a TOML run file')
    try:
        runs = [run_config.load(path, SCHEMA) for path in args.runs] if args.runs else [
            run_config.resolve(SCHEMA, parameters={'BETA': beta/10, 'L': length},
                               output={'results': f'ising.L_{length}beta_{beta/10}.h5'})
            for beta in range(11) for length in sizes]
        inputs = {Path(path).resolve() for path in args.runs}
        outputs = [Path(run.output['results']).resolve() for run in runs]
        if len(set(outputs)) != len(outputs) or inputs.intersection(outputs):
            raise ValueError('Output paths must be distinct and must not replace run files')
        simulations = [simulation(run.parameters['BETA'], run.parameters['L'],
                                  run.execution['seed'], run.execution['bins']) for run in runs]
        for run, output, sim in zip(runs, outputs, simulations):
            if args.validate:
                print(f'Valid Ising configuration: {output}')
                continue
            sim.run(run.parameters['THERMALIZATION'], run.parameters['SWEEPS'])
            sim.save(output)
        if plot and not args.validate:
            import matplotlib.pyplot as plt
            import pyalps.plot
            data = pyalps.loadMeasurements([str(path) for path in outputs], ['Binder Ratio'])
            pyalps.plot.plot(pyalps.collectXY(data, x='BETA', y='Binder Ratio', foreach=['L']))
            plt.xlabel(r'Inverse Temperature $\beta$')
            plt.ylabel(r'Binder ratio $\langle m^4\rangle/\langle m^2\rangle^2$')
            plt.title('2D Ising model')
            plt.legend()
            plt.show()
    except (ValueError, RuntimeError, OSError) as error:
        parser.exit(1, f'{parser.prog}: {error}\n')


if __name__ == '__main__':
    main()
