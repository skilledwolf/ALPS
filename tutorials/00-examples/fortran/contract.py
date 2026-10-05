# Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
"""Installed Fortran bridge: independent chains, physics and exact restart."""
from pathlib import Path
import itertools
import shutil
import subprocess
import sys
import tempfile
import tomllib
import unittest

import h5py
import numpy as np


class FortranContract(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="alps-fortran-")
        self.addCleanup(temporary.cleanup)
        self.directory = Path(temporary.name)

    def run_case(self, name, *, seed=42, chains=2, sweeps=257, warmup=13,
                 budget=0, checkpoint=None, parameters="", validate=False, success=True):
        path = self.directory / (name + ".toml")
        if self.kind == "ising":
            parameters = f"L=2\nTEMPERATURE=3.0\nSWEEPS={sweeps}\nTHERMALIZATION={warmup}\n" + parameters
        text = f'''[parameters]
{parameters}
[input]
{f'checkpoint = "{checkpoint}"' if checkpoint else ''}
[output]
results = "{name}.h5"
checkpoint = "{name}.checkpoint.h5"
[execution]
seed = {seed}
chains = {chains}
bins = 8
max_sweeps = {budget}
'''
        path.write_text(text)
        result = subprocess.run([str(self.executable), *( ["--validate"] if validate else []), str(path)],
                                capture_output=True, text=True, timeout=60)
        self.assertEqual(result.returncode == 0, success, result.stdout + result.stderr)
        return path.with_suffix(".checkpoint.h5")

    @staticmethod
    def datasets(path, prefix):
        values = {}
        with h5py.File(path) as archive:
            archive[prefix].visititems(lambda name, value: values.update({name: value[()]})
                                      if isinstance(value, h5py.Dataset) else None)
        return values

    def same_state(self, first, second, clone1=0, clone2=0):
        prefix = "simulation/realizations/0/clones/"
        a = self.datasets(first, prefix + str(clone1))
        b = self.datasets(second, prefix + str(clone2))
        # Identity differs when comparing a pooled run against separate runs.
        a.pop("checkpoint/fortran/chain"); b.pop("checkpoint/fortran/chain")
        self.assertEqual(a.keys(), b.keys())
        for name in a:
            np.testing.assert_array_equal(a[name], b[name], err_msg=name)

    def test_schema_and_validation(self):
        result = subprocess.run([str(self.executable), "--schema"], capture_output=True, text=True, check=True)
        self.assertIn("checkpoint", tomllib.loads(result.stdout)["input"])
        self.run_case("validate", validate=True)
        self.assertFalse(list(self.directory.glob("*.h5")))

    def test_independent_chains(self):
        multiple = self.run_case("multiple")
        for chain in range(2):
            separate = self.run_case(f"separate{chain}", chains=1, seed=42+chain)
            self.same_state(multiple, separate, chain, 0)
        if self.kind == "ising":
            with h5py.File(multiple) as ar:
                root = ar["simulation/realizations/0/clones"]
                self.assertNotEqual(root['0/checkpoint/engine/engine'][()], root['1/checkpoint/engine/engine'][()])

    def test_exact_restart(self):
        full = self.run_case("full")
        for budget in (5, 71):
            with self.subTest(budget=budget):
                part = self.run_case(f"part{budget}", budget=budget)
                resumed = self.run_case(f"resumed{budget}", checkpoint=part)
                for chain in range(2):
                    self.same_state(full, resumed, chain, chain)

    def test_extend_production(self):
        short = self.run_case("short", sweeps=31)
        extended = self.run_case("extended", checkpoint=short)
        full = self.run_case("full")
        for chain in range(2):
            self.same_state(full, extended, chain, chain)

    def test_reject_corruption_without_outputs(self):
        source = self.run_case("source", budget=5)
        prefix = "simulation/realizations/0/clones/0/checkpoint/fortran/"
        for kind in ("type", "count", "fields", "version", "state"):
            with self.subTest(kind=kind):
                corrupt = self.directory / f"corrupt-{kind}.h5"
                shutil.copyfile(source, corrupt)
                with h5py.File(corrupt, "r+") as ar:
                    if kind in ("type", "count"):
                        ar[prefix+"fields/0"].attrs[kind] = 12345
                    elif kind == "fields":
                        ar[prefix+"fields_count"][()] = 12345
                    elif kind == "version":
                        ar[prefix+"version"][()] = 12345
                    else:
                        del ar[prefix+"fields/1/value"]
                self.run_case("bad", checkpoint=corrupt, success=False)
                self.assertFalse((self.directory/"bad.h5").exists())
                self.assertFalse((self.directory/"bad.checkpoint.h5").exists())

    def test_physics_or_character_error(self):
        if self.kind == "hello":
            self.run_case("long", parameters='WORLD="' + 'x'*100 + '"', success=False)
            self.run_case("values", parameters='X=2.5')
            with h5py.File(self.directory/"values.h5") as ar:
                np.testing.assert_array_equal(ar["simulation/results/Values/mean/value"][()], [2.5, 5., 7.5, 10.])
            return
        self.run_case("physics", sweeps=100000, warmup=1000)
        energies = []
        for values in itertools.product((-1, 1), repeat=4):
            spins = np.array(values).reshape(2, 2)
            energies.append(-np.sum(spins * (np.roll(spins, 1, 0) + np.roll(spins, 1, 1))))
        energies = np.array(energies)
        weights = np.exp(-energies/3.0)
        exact = np.sum(weights*energies)/np.sum(weights)/4
        with h5py.File(self.directory/"physics.h5") as ar:
            result = ar["simulation/results/Energy"]
            self.assertAlmostEqual(float(np.asarray(result['mean/value'][()]).item()), exact, delta=0.025)


if __name__ == "__main__":
    FortranContract.executable = Path(sys.argv.pop(1)).resolve()
    FortranContract.kind = sys.argv.pop(1)
    unittest.main()
