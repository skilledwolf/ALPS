"""Check TOML orchestration, validation, and fail-before-output behavior."""
import argparse
from pathlib import Path
import subprocess
import tempfile
import tomllib
import unittest

class CTIntCLIContract(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="alps-ctint-cli-")
        self.addCleanup(temporary.cleanup)
        self.directory = Path(temporary.name)

    def invoke(self, *args, unchanged=True):
        before = {path.name: path.read_bytes() for path in self.directory.iterdir() if path.is_file()}
        result = subprocess.run([str(self.executable), *map(str, args)], cwd=self.directory,
                                capture_output=True, text=True, timeout=15)
        if unchanged:
            self.assertEqual({path.name: path.read_bytes() for path in self.directory.iterdir() if path.is_file()}, before)
        return result

    def config(self):
        path = self.directory / "run.toml"
        path.write_text('''[parameters]
BETA=2.0
U=0.0
MU=0.0
ALPHA=-0.01
N=8
NMATSUBARA=4
SWEEPS=8
THERMALIZATION=2
MEASUREMENT_PERIOD=1
[input]
atomic=true
[output]
results="result.h5"
[execution]
time_limit=0
''')
        return path

    def test_help_and_invalid_cli(self):
        result = self.invoke("--help", "missing.toml")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("--validate", result.stdout)
        for args in [(), ("--bad-option",), ("missing.toml",), ("a.toml", "b.toml")]:
            with self.subTest(args=args):
                self.assertNotEqual(self.invoke(*args).returncode, 0)

    def test_validation_has_no_output(self):
        result = self.invoke("--validate", self.config())
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_multiband_schema_and_validation(self):
        path = self.config()
        path.write_text(path.read_text().replace("N=8", "N=8\nFLAVORS=4\nEPS_3=0.0\nEPSSQ_3=0.0"))
        result = self.invoke("--schema", path)
        self.assertEqual(result.returncode, 0, result.stderr)
        schema = tomllib.loads(result.stdout)
        self.assertIn("EPS_3", schema["parameters"])
        self.assertNotIn("EPS_4", schema["parameters"])
        self.assertIn("interaction_matrix", schema["input"])
        result = self.invoke("--validate", path)
        self.assertEqual(result.returncode, 0, result.stderr)
        path.write_text(path.read_text().replace("EPS_3=0.0", "EPS_4=0.0"))
        self.assertNotEqual(self.invoke("--validate", path).returncode, 0)
        for value in ("true", "3.0", "1", "129"):
            path = self.config()
            path.write_text(path.read_text().replace("N=8", f"N=8\nFLAVORS={value}"))
            self.assertNotEqual(self.invoke("--schema", path).returncode, 0)

    def test_explicit_matrix_and_sparse_rows(self):
        path = self.config()
        matrix = self.directory / "interaction.dat"
        path.write_text(path.read_text().replace("\nU=0.0\n", "\nFLAVORS=3\n")
                        .replace("atomic=true", 'atomic=true\ninteraction_matrix="interaction.dat"'))
        for text in ("", "0 1 1.0\n1 0 1.0\n"):
            matrix.write_text(text)
            result = self.invoke("--validate", path)
            self.assertEqual(result.returncode, 0, result.stderr)
        for text, message in (("0 1 1.0\n", "symmetric"), ("0 0 1.0\n", "zero diagonal"),
                              ("3 0 1.0\n", "Invalid index"), ("0 1 broken\n", "Malformed")):
            matrix.write_text(text)
            result = self.invoke("--validate", path)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn(message, result.stderr)

    def test_invalid_settings_preserve_outputs(self):
        changes = [("N=8", "N=0"), ("N=8", "N_TAU=8"),
                   ("atomic=true", "atomic=false"),
                   ("time_limit=0", "time_limit=-1"),
                   ("NMATSUBARA=4", "NMATSUBARA=4\nNMATSUBARA_MEASUREMENTS=5")]
        for old, new in changes:
            with self.subTest(new=new):
                path = self.config()
                path.write_text(path.read_text().replace(old, new))
                (self.directory / "result.h5").write_bytes(b"existing scientific output")
                self.assertNotEqual(self.invoke(path).returncode, 0)

    def test_output_cannot_replace_run_file(self):
        path = self.config()
        path.write_text(path.read_text().replace('results="result.h5"', 'results="run.toml"'))
        self.assertNotEqual(self.invoke(path).returncode, 0)

    def test_legacy_parameter_input_rejected(self):
        path = self.directory / "legacy.parm"
        path.write_text("BETA=2; U=0; N=8;")
        self.assertNotEqual(self.invoke(path).returncode, 0)

    def test_small_atomic_run_finishes_and_writes_only_results(self):
        result = self.invoke(self.config(), unchanged=False)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual((self.directory / "result.h5").read_bytes()[:8], b"\x89HDF\r\n\x1a\n")
        self.assertEqual({path.name for path in self.directory.iterdir()}, {"run.toml", "result.h5"})

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("interaction", type=Path)
    arguments = parser.parse_args()
    CTIntCLIContract.executable = arguments.interaction.resolve(strict=True)
    result = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(CTIntCLIContract))
    raise SystemExit(0 if result.wasSuccessful() else 1)
