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

from contextlib import nullcontext
import numpy as np
import pyalps.hdf5 as hdf5
import pyalps.ngs as ngs
from pyalps.alea import BatchAccumulator
from pyalps import hdf5_name_encode


def configuration(parameters):
    length, thermalization, total = (parameters[key] for key in ('L', 'THERMALIZATION', 'SWEEPS'))
    temperature = parameters['T']
    if (any(not isinstance(value, (int, np.integer)) or isinstance(value, (bool, np.bool_))
            for value in (length, thermalization, total))
            or length < 2 or thermalization < 0 or total <= 0
            or not isinstance(temperature, (int, float, np.integer, np.floating))
            or isinstance(temperature, (bool, np.bool_))
            or not np.isfinite(temperature) or temperature <= 0):
        raise ValueError('invalid Ising configuration')
    beta = 1. / float(temperature)
    if not np.isfinite(beta):
        raise ValueError('invalid Ising inverse temperature')
    return int(length), int(thermalization), int(total), beta


class sim(ngs.mcbase):
    def __init__(self, parameters, seed=42):
        super().__init__(parameters, seed)
        self.length, self.thermalization_sweeps, self.total_sweeps, self.beta = configuration(self.parameters)
        self.sweeps = 0
        self.spins = np.array([-1 if self.random() < 0.5 else 1 for _ in range(self.length)])
        for name in ('Energy', 'Magnetization', 'Magnetization^2', 'Magnetization^4'):
            self.measurements[name] = BatchAccumulator(size=1, num_batches=64)
        self.measurements['Correlations'] = BatchAccumulator(size=self.length, num_batches=64)

    def update(self):
        for _ in range(self.length):
            i = int(self.length * self.random())
            neighbors = self.spins[(i + 1) % self.length] + self.spins[(i - 1) % self.length]
            probability = np.exp(-2 * self.beta * self.spins[i] * neighbors)
            if probability >= 1 or self.random() < probability:
                self.spins[i] = -self.spins[i]

    def measure(self):
        self.sweeps += 1
        if self.sweeps > self.thermalization_sweeps:
            energy = -np.dot(self.spins, np.roll(self.spins, 1)) / self.length
            magnetization = self.spins.mean()
            correlations = [np.dot(self.spins, np.roll(self.spins, distance)) / self.length
                            for distance in range(self.length)]
            self.measurements['Energy'] << energy
            self.measurements['Magnetization'] << magnetization
            self.measurements['Magnetization^2'] << magnetization**2
            self.measurements['Magnetization^4'] << magnetization**4
            self.measurements['Correlations'] << correlations

    def fraction_completed(self):
        return max(self.sweeps - self.thermalization_sweeps, 0) / self.total_sweeps

    def save(self, ar):
        with ar.native() if isinstance(ar, hdf5.archive) else nullcontext(ar) as native:
            ngs.mcbase.save(self, native)
            native['checkpoint/sweeps'] = self.sweeps
            native['checkpoint/spins'] = self.spins

    def load(self, ar):
        with ar.native() if isinstance(ar, hdf5.archive) else nullcontext(ar) as native:
            parameters = ngs.params(native, '/parameters')
            length, thermalization, total, beta = configuration(parameters)
            sweeps = native['checkpoint/sweeps']
            spins = native['checkpoint/spins']
            if (not isinstance(sweeps, np.integer) or not 0 <= sweeps <= thermalization + total
                    or spins.shape != (length,) or spins.dtype.kind not in 'if'
                    or not np.all((spins == -1) | (spins == 1))):
                raise ValueError('invalid Ising checkpoint')
            count = max(int(sweeps) - thermalization, 0)
            for name in self.measurements:
                counts = native['measurements/' + hdf5_name_encode(name) + '/batch/count']
                if sum(int(value) for value in counts) != count:
                    raise ValueError('Ising measurement count does not match checkpoint progress')
            ngs.mcbase.load(self, native)
        self.length = length
        self.sweeps = int(sweeps)
        self.thermalization_sweeps = thermalization
        self.total_sweeps = total
        self.beta = beta
        self.spins = spins
