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

import pyalps
import pyalps.hdf5 as hdf5
from pyalps.alea import BatchAccumulator
from pyalps.ngs import params, random01
from contextlib import nullcontext
import numpy as np


class sim:

    def __init__(self, parameters):
        self.parameters = params(parameters)
        self.random = random01(self.parameters.valueOrDefault('SEED', 42))
        self.length = int(self.parameters['L'])
        self.sweeps = 0
        self.thermalization_sweeps = int(self.parameters['THERMALIZATION'])
        self.total_sweeps = int(self.parameters['SWEEPS'])
        self.beta = 1. / float(self.parameters['T'])
        self.spins = np.array([(-x if self.random() < 0.5 else x) for x in np.ones(self.length)])
        self.measurements = {
            name: BatchAccumulator() for name in
            ('Energy', 'Magnetization', 'Magnetization^2', 'Magnetization^4')
        }
        self.measurements['Correlations'] = BatchAccumulator(self.length)

    def update(self):
        for j in range(self.length):
            i = int(float(self.length) * self.random())
            right = i + 1 if i + 1 < self.length else 0
            left = self.length - 1 if i - 1 < 0 else i - 1
            p = np.exp(-2. * self.beta * self.spins[i] * (self.spins[right] + self.spins[left]))
            if p >= 1. or self.random() < p:
                self.spins[i] =- self.spins[i]

    def measure(self):
        self.sweeps += 1
        if self.sweeps > self.thermalization_sweeps:
            tmag = 0
            ten = 0
            corr = np.zeros(self.length)
            for i in range(self.length):
                tmag += self.spins[i]
                ten += -self.spins[i] * self.spins[i + 1 if i + 1 < self.length else 0]
            for d in range(self.length):
                corr[d] = np.inner(self.spins, np.roll(self.spins, d)) / float(self.length)
            ten /= self.length
            tmag /= self.length
            self.measurements['Energy'] << ten
            self.measurements['Magnetization'] << tmag
            self.measurements['Magnetization^2'] << tmag**2
            self.measurements['Magnetization^4'] << tmag**4
            self.measurements['Correlations'] << corr

    def fraction_completed(self):
        return 0 if self.sweeps < self.thermalization_sweeps else (self.sweeps - self.thermalization_sweeps) / float(self.total_sweeps)

    def run(self, stopCallback):
        while self.fraction_completed() < 1.:
            self.update()
            self.measure()
            if stopCallback():
                return False
        return True

    def collectResults(self):
        return {name: accumulator.result() for name, accumulator in self.measurements.items()}

    def save(self, ar):
        # Atomic checkpoint callbacks already own a NativeArchive. Ordinary
        # Python archives transfer ownership once for the complete operation.
        with ar.native() if isinstance(ar, hdf5.archive) else nullcontext(ar) as native:
            native['/parameters'] = self.parameters
            base = '/simulation/realizations/0/clones/0'
            for name, accumulator in self.measurements.items():
                accumulator.save(native, base + '/measurements/' + pyalps.hdf5_name_encode(name))
            native[base + '/checkpoint/sweeps'] = self.sweeps
            native[base + '/checkpoint/spins'] = self.spins
            native[base + '/checkpoint'] = self.random

    def load(self, ar):
        with ar.native() if isinstance(ar, hdf5.archive) else nullcontext(ar) as native:
            restored = sim(dict(params(native, '/parameters')))
            base = '/simulation/realizations/0/clones/0'
            for name, accumulator in restored.measurements.items():
                size = accumulator.size
                accumulator.load(native, base + '/measurements/' + pyalps.hdf5_name_encode(name))
                if accumulator.size != size:
                    raise ValueError('invalid Ising measurement shape')
            sweeps = native[base + '/checkpoint/sweeps']
            restored.spins = native[base + '/checkpoint/spins']
            if (not isinstance(sweeps, np.integer) or sweeps < 0
                    or sweeps > restored.thermalization_sweeps + restored.total_sweeps
                    or restored.spins.shape != (restored.length,)
                    or not np.all(np.abs(restored.spins) == 1)):
                raise ValueError('invalid Ising checkpoint')
            restored.sweeps = int(sweeps)
            count = max(restored.sweeps - restored.thermalization_sweeps, 0)
            if any(accumulator.count != count for accumulator in restored.measurements.values()):
                raise ValueError('Ising measurement count does not match checkpoint progress')
            context = native.context
            try:
                native.set_context(base + '/checkpoint')
                restored.random.load(native)
            finally:
                native.set_context(context)
        self.__dict__ = restored.__dict__
