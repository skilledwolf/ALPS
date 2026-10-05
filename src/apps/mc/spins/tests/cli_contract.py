# Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
"""Native validation must finish before replacing any scientific output."""
from pathlib import Path
import subprocess
import sys
import tempfile
import tomllib
import unittest


class SpinMCContract(unittest.TestCase):
    def test_schema_and_validation_rollback(self):
        with tempfile.TemporaryDirectory(prefix="alps-spinmc-cli-") as temporary:
            directory = Path(temporary)
            path = directory / "run.toml"
            output = directory / "result.h5"
            output.write_bytes(b"existing scientific analysis")
            valid = '''[parameters]
MODEL = "XY"
UPDATE = "local"
LATTICE = "chain lattice"
L = 4
T = 1.6
J = [1.0, 0.2, -0.1, 0.8]
J0 = [0.9, 0.7]
D0 = [0.1, 0.2]
S0 = 1.5
h = [0.1, -0.2]
SWEEPS = 37
[execution]
chains = 2
bins = 8
seed = 2873
[output]
results = "result.h5"
'''
            path.write_text(valid)
            schema = subprocess.run([self.executable, "--schema", str(path)], capture_output=True, text=True)
            self.assertEqual(schema.returncode, 0, schema.stderr)
            definitions = tomllib.loads(schema.stdout)
            self.assertEqual(definitions["parameters"]["J0"]["type"], "float64[]")
            self.assertEqual(definitions["parameters"]["S0"]["type"], "float64")
            self.assertEqual(definitions["execution"]["print_sweeps"]["type"], "int64")
            changes = [(None, None), ('T = 1.6', 'T = 0.0'),
                       ('T = 1.6', 'T = 1.6\nbeta = 0.5'),
                       ('UPDATE = "local"', 'UPDATE = "cluster"'),
                       ('J0 = [0.9, 0.7]', 'J0 = "0.9 0.7"'),
                       ('L = 4', 'L = 4\nDISORDERSEED = 3'),
                       ('bins = 8', 'bins = 3'),
                       ('seed = 2873', 'seed = 2147483647'),
                       ('seed = 2873', 'seed = 2873\nrng = "unknown_rng"'),
                       ('seed = 2873', 'seed = 2873\nerror_variable = "Energy"'),
                       ('results = "result.h5"', 'results = "run.toml"')]
            for old, new in changes:
                with self.subTest(change=new):
                    path.write_text(valid if old is None else valid.replace(old, new))
                    before = {p.name: p.read_bytes() for p in directory.iterdir()}
                    result = subprocess.run([self.executable, "--validate", str(path)],
                                            capture_output=True, text=True, timeout=20)
                    self.assertEqual(result.returncode == 0, old is None, result.stderr)
                    self.assertEqual({p.name: p.read_bytes() for p in directory.iterdir()}, before)


if __name__ == "__main__":
    SpinMCContract.executable = str(Path(sys.argv[1]).resolve(strict=True))
    result = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(SpinMCContract))
    raise SystemExit(not result.wasSuccessful())
