# Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
"""Exercise the TOML-only Hirsch-Fye CLI without modifying scientific inputs."""
import argparse
from pathlib import Path
import subprocess
import tempfile
import unittest


class HirschFyeCLIContract(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="alps-hirschfye-cli-")
        self.addCleanup(temporary.cleanup)
        self.directory = Path(temporary.name)

    def invoke(self, *arguments, unchanged=True):
        before = {path.relative_to(self.directory): path.read_bytes()
                  for path in self.directory.rglob("*") if path.is_file()}
        result = subprocess.run([str(self.executable), *map(str, arguments)],
                                cwd=self.directory, capture_output=True,
                                text=True, timeout=15)
        if unchanged:
            self.assertEqual({path.relative_to(self.directory): path.read_bytes()
                              for path in self.directory.rglob("*") if path.is_file()}, before)
        return result

    def config(self):
        path = self.directory / "run.toml"
        path.write_text('''[parameters]
BETA=2.0
U=0.0
EPSSQ_0=0.0
EPSSQ_1=0.0
N=4
NMATSUBARA=4
SWEEPS=4
THERMALIZATION=1
[input]
g0="g0.h5"
[output]
results="result.h5"
[execution]
time_limit=0
''')
        return path

    def fixture(self, kind):
        result = subprocess.run([str(self.fixture_writer), str(self.directory / "g0.h5"), kind],
                                capture_output=True, text=True, timeout=15)
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_valid_configuration_checks_input_without_output(self):
        self.fixture("valid")
        (self.directory / "result.h5").write_bytes(b"existing scientific results")
        result = self.invoke("--validate", self.config())
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("Valid Hirsch-Fye", result.stdout)

    def test_malformed_green_functions_preserve_output(self):
        for kind in ["short", "oversized", "real", "missing-flavor", "nonfinite"]:
            with self.subTest(kind=kind):
                self.fixture(kind)
                path = self.config()
                (self.directory / "result.h5").write_bytes(b"existing scientific results")
                for arguments in [(path,), ("--validate", path)]:
                    self.assertNotEqual(self.invoke(*arguments).returncode, 0)

    def test_free_solver_preserves_input_and_matches_exact_result(self):
        self.fixture("valid")
        path = self.config()
        source = (self.directory / "g0.h5").read_bytes()
        document = path.read_bytes()
        result = self.invoke(path, unchanged=False)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual((self.directory / "g0.h5").read_bytes(), source)
        self.assertEqual(path.read_bytes(), document)
        self.assertEqual({item.name for item in self.directory.iterdir()},
                         {"g0.h5", "run.toml", "result.h5"})
        checked = subprocess.run([str(self.fixture_writer), str(self.directory / "result.h5"),
                                  "verify-output"], capture_output=True, text=True, timeout=15)
        self.assertEqual(checked.returncode, 0, checked.stderr)

    def test_help_and_cli_errors(self):
        result = self.invoke("--help", "missing.toml")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("--validate", result.stdout)
        for arguments in [(), ("--bad-option",), ("missing.toml",), ("a.toml", "b.toml")]:
            with self.subTest(arguments=arguments):
                result = self.invoke(*arguments)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("hirschfye:", result.stderr)

    def test_missing_or_invalid_g0_preserves_output(self):
        path = self.config()
        (self.directory / "result.h5").write_bytes(b"existing scientific results")
        for content in [None, b"not an HDF5 archive"]:
            with self.subTest(content=content):
                if content is not None:
                    (self.directory / "g0.h5").write_bytes(content)
                for arguments in [(path,), ("--validate", path)]:
                    self.assertNotEqual(self.invoke(*arguments).returncode, 0)

    def test_invalid_settings_preserve_output(self):
        changes = [("BETA=2.0", "BETA=0.0"), ("U=0.0", "U=-1.0"),
                   ("N=4", "N=4.5"), ("NMATSUBARA=4", "NMATSUBARA=0"),
                   ("SWEEPS=4", "SWEEPS=-1"), ("U=0.0", "U=0.0\nFLAVORS=4"),
                   ("time_limit=0", 'time_limit=0\nloop="tau"'),
                   ("time_limit=0", 'time_limit=0\nsolver="interaction"'),
                   ("time_limit=0", "time_limit=-1"),
                   ('g0="g0.h5"', "")]
        self.fixture("valid")
        (self.directory / "result.h5").write_bytes(b"existing scientific results")
        for old, new in changes:
            with self.subTest(new=new):
                path = self.config()
                path.write_text(path.read_text().replace(old, new))
                self.assertNotEqual(self.invoke("--validate", path).returncode, 0)

    def test_output_cannot_replace_run_or_input(self):
        self.fixture("valid")
        for name in ["run.toml", "g0.h5"]:
            with self.subTest(name=name):
                path = self.config()
                path.write_text(path.read_text().replace('results="result.h5"',
                                                        f'results="{name}"'))
                self.assertNotEqual(self.invoke(path).returncode, 0)

    def test_schema_is_the_installed_contract(self):
        result = self.invoke("--schema")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('application = "hirschfye"', result.stdout)

    def test_legacy_run_formats_are_rejected(self):
        self.fixture("valid")
        for name, content in [("legacy.parm", b"BETA=2; U=0; N=8;"),
                              ("legacy.h5", (self.directory / "g0.h5").read_bytes())]:
            with self.subTest(name=name):
                path = self.directory / name
                path.write_bytes(content)
                self.assertNotEqual(self.invoke(path).returncode, 0)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("hirschfye", type=Path)
    parser.add_argument("fixture_writer", type=Path)
    arguments = parser.parse_args()
    HirschFyeCLIContract.executable = arguments.hirschfye.resolve(strict=True)
    HirschFyeCLIContract.fixture_writer = arguments.fixture_writer.resolve(strict=True)
    result = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(HirschFyeCLIContract))
    raise SystemExit(0 if result.wasSuccessful() else 1)
