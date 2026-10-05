# Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
"""Validate the native simplemc entry point without a Python binding build."""
from pathlib import Path
import subprocess
import sys
import tempfile
import tomllib
import unittest


class SimpleMCContract(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="alps-simplemc-cli-")
        self.addCleanup(temporary.cleanup)
        self.directory = Path(temporary.name)

    def config(self, model="ising", name=None):
        path = self.directory / (name or model + ".toml")
        path.write_text(Path(__file__).with_name(model + ".toml").read_text(), encoding="utf-8")
        return path

    def invoke(self, *args, unchanged=True):
        before = {p.name: p.read_bytes() for p in self.directory.iterdir() if p.is_file()}
        result = subprocess.run([str(self.executable), *map(str, args)], cwd=self.directory,
                                text=True, capture_output=True, timeout=20)
        if unchanged:
            self.assertEqual({p.name: p.read_bytes() for p in self.directory.iterdir() if p.is_file()}, before)
        return result

    def test_validation_and_dynamic_graph_schema(self):
        for model in ("ising", "xy", "heisenberg"):
            with self.subTest(model=model):
                path = self.config(model)
                result = self.invoke("--validate", path)
                self.assertEqual(result.returncode, 0, result.stderr)
                result = self.invoke("--schema", path)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(tomllib.loads(result.stdout)["parameters"]["L"]["type"], "int64")

    def test_invalid_settings_preserve_outputs(self):
        changes = [('T = 2.0', 'T = 0.0'), ('chains = 2', 'chains = 0'),
                   ('bins = 8', 'bins = 3'), ('seed = 2873', 'seed = 2147483647'),
                   ('seed = 2873', 'seed = 2873\nrng = "unknown_rng"'),
                   ('L = 4', 'L = 4\nNUM_CLONES = 2'),
                   ('L = 4', 'L = 4\nDISORDERSEED = 3'),
                   ('SWEEPS = 37', 'SWEEPS = "[37:]"'),
                   ('results = "ising.h5"', 'results = "ising.toml"')]
        for old, new in changes:
            with self.subTest(new=new):
                path = self.config()
                path.write_text(path.read_text().replace(old, new), encoding="utf-8")
                (self.directory / "ising.h5").write_bytes(b"existing scientific data")
                result = self.invoke(path)
                self.assertNotEqual(result.returncode, 0)

    def test_implicit_lattice_library_is_protected(self):
        library = self.directory / "lattices.xml"
        library.write_bytes((Path(__file__).resolve().parents[4] /
                             "alps/resources/lattices.xml").read_bytes())
        for key in ("results", "checkpoint"):
            path = self.config()
            text = path.read_text()
            if key == "results":
                text = text.replace('results = "ising.h5"', 'results = "lattices.xml"')
            else:
                text = text.replace('[output]', '[output]\ncheckpoint = "lattices.xml"')
            path.write_text(text)
            result = self.invoke(path)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("replace", result.stderr)

    def test_multiple_tasks_preflight_and_snapshot_namespace(self):
        first, second = self.config(), self.config("xy")
        second.write_text(second.read_text().replace("T = 1.6", "T = -1.0"))
        self.assertNotEqual(self.invoke(first, second).returncode, 0)
        second = self.config("xy")
        first.write_text(first.read_text().replace('results = "ising.h5"',
                        'results = "ising.h5"\nsnapshot_prefix = "snap"'))
        second.write_text(second.read_text().replace('results = "xy.h5"',
                          'results = "snap.clone1.8.vtk"'))
        result = self.invoke(first, second)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Snapshot prefix", result.stderr)

    def test_multiple_models_write_native_hdf5(self):
        files = [self.config(model) for model in ("ising", "xy", "heisenberg")]
        result = self.invoke(*files, unchanged=False)
        self.assertEqual(result.returncode, 0, result.stderr)
        for model in ("ising", "xy", "heisenberg"):
            self.assertEqual((self.directory / (model + ".h5")).read_bytes()[:8], b"\x89HDF\r\n\x1a\n")

    def test_help_invalid_cli_and_legacy_input(self):
        self.assertEqual(self.invoke("--help").returncode, 0)
        for args in ((), ("--unknown",), ("missing.toml",)):
            self.assertNotEqual(self.invoke(*args).returncode, 0)
        legacy = self.directory / "old.parm"
        legacy.write_text('LATTICE="chain lattice"; L=4; SWEEPS=37;')
        self.assertNotEqual(self.invoke(legacy).returncode, 0)


if __name__ == "__main__":
    SimpleMCContract.executable = Path(sys.argv[1]).resolve(strict=True)
    result = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(SimpleMCContract))
    raise SystemExit(not result.wasSuccessful())
