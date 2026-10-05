# Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
"""Real DMFT -> CT-QMC -> scientific results contracts; Python stdlib only."""
import argparse
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile
import unittest


class DMFTCLIContract(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="alps-dmft-cli-")
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name) / "run files"
        self.directory.mkdir()
        self.work = Path(self.temporary.name) / "unrelated working directory"
        self.work.mkdir()

    def invoke(self, *arguments):
        # Built-in solvers are found in the ALPS bin directory.
        environment = dict(os.environ, ALPS_BIN_PATH=str(self.hybridization.parent))
        return subprocess.run([str(self.dmft), *map(str, arguments)], cwd=self.work, env=environment,
                              capture_output=True, text=True, timeout=90)

    def config(self, solver="interaction", loop="omega", text=False, retarded=False):
        if text:
            (self.directory / "text").mkdir(exist_ok=True)
        specific, sweeps, thermalization = {
            "hybridization": ("N_MEAS=4\nMEASURE_freq=true\nMEASURE_time=true\n",
                              3000, 100),
            "hirschfye": ("", 200, 20),
        }.get(solver, ("ALPHA=-0.01\nMEASUREMENT_PERIOD=1\n", 32, 4))
        retarded_settings = ""
        if retarded:
            (self.directory / "retarded.dat").write_text("".join(f"{i / 8:.17g} 0 0\n" for i in range(17)))
            retarded_settings = ('retarded_interaction="retarded.dat"\n'
                                 'retarded_interaction_format="text"\n'
                                 'retarded_interaction_coordinate="tau"\n')
        path = self.directory / "run.toml"
        path.write_text(f'''[parameters]
BETA=2.0
U=0.0
MU=0.0
H=0.0
t=1.0
N=16
NMATSUBARA=8
FLAVORS=2
SITES=1
SWEEPS={sweeps}
THERMALIZATION={thermalization}
CONVERGED=0.0
SYMMETRIZATION=true
{specific}[input]
{retarded_settings}[output]
results="results.h5"
text={str(text).lower()}
text_directory="text"
[execution]
solver={json.dumps(solver)}
loop={json.dumps(loop)}
max_iterations=2
seed=42
time_limit=0
''', encoding="utf-8")
        return path

    def verify(self, mode):
        """mode is "exact" (U=0 equals G0), "sampled" (analytic within noise) or "retarded"."""
        result = subprocess.run([str(self.helper), "--verify", str(self.directory / "results.h5"), mode],
                                cwd=self.work, capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def run_case(self, solver, reference, loop="omega", text=False, retarded=False):
        path = self.config(solver, loop, text, retarded)
        original = path.read_bytes()
        original_retarded = (self.directory / "retarded.dat").read_bytes() if retarded else None
        validated = self.invoke("--validate", path)
        self.assertEqual(validated.returncode, 0, validated.stdout + validated.stderr)
        self.assertFalse((self.directory / "results.h5").exists())
        self.assertEqual(path.read_bytes(), original)
        self.assertEqual(list(self.work.iterdir()), [])
        result = self.invoke(path)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(path.read_bytes(), original)
        self.assertEqual(list(self.work.iterdir()), [], "DMFT wrote implicit working-directory files")
        self.assertTrue((self.directory / "results.h5").is_file(), "relative output did not resolve from run file")
        expected = {"run.toml", "results.h5"}
        if retarded:
            expected.add("retarded.dat")
            self.assertEqual((self.directory / "retarded.dat").read_bytes(), original_retarded)
        if text:
            expected.add("text")
            written = [item.name for item in (self.directory / "text").iterdir()]
            self.assertTrue(written, "requested text output was not produced")
            # Driver text settings are not forwarded to the impurity solver.
            self.assertFalse([name for name in written if name.endswith(".dat")], written)
        self.assertEqual({p.name for p in self.directory.iterdir()}, expected)
        self.verify(reference)

    def test_cthyb_two_omega_iterations(self):
        self.run_case("hybridization", "sampled")

    def test_cthyb_two_time_iterations_with_text_output(self):
        self.run_case("hybridization", "sampled", loop="tau", text=True)

    def test_cthyb_zero_retarded_kernel_retains_the_u0_reference(self):
        self.run_case("hybridization", "sampled", retarded=True)

    def test_hirsch_fye_two_iterations(self):
        self.run_case("hirschfye", "sampled")

    def test_hirsch_fye_validation_uses_the_standalone_numerical_domain(self):
        changes = [("U=0.0", "U=1.0e300"),
                   ("U=0.0", "U=0.0\nMEASURE_FOURPOINT_FUNCTION=false"),
                   ("seed=42", "seed=42\nbins=3"),
                   ('solver="hirschfye"', 'solver="Hirsch-Fye"')]
        for old, new in changes:
            with self.subTest(new=new):
                path = self.config("hirschfye")
                path.write_text(path.read_text().replace(old, new))
                self.validate_preserving_files(path)

    def test_multiband_ctint_forwards_the_interaction_matrix(self):
        path = self.config()
        path.write_text(path.read_text().replace("FLAVORS=2", "FLAVORS=4")
                        .replace("\nU=0.0\n", "\nU=0.5\n")
                        .replace("[input]", '[input]\ninteraction_matrix="zero matrix.dat"'))
        matrix = self.directory / "zero matrix.dat"
        matrix.write_text("0 1 0\n1 0 0\n")
        before = path.read_bytes(), matrix.read_bytes()
        validated = self.invoke("--validate", path)
        self.assertEqual(validated.returncode, 0, validated.stdout + validated.stderr)
        result = self.invoke(path)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(before, (path.read_bytes(), matrix.read_bytes()))
        self.assertEqual(list(self.work.iterdir()), [])
        self.assertEqual({p.name for p in self.directory.iterdir()},
                         {"run.toml", "zero matrix.dat", "results.h5"})
        # Scalar U is nonzero: dropping the zero matrix changes this exact result.
        self.verify("exact")

    def test_removed_scheduler_selector_requires_a_custom_solver_contract(self):
        path = self.config("Interaction Expansion")
        self.validate_preserving_files(path)

    def test_multiband_matrix_validation_matches_the_standalone_solver(self):
        for content in ("0 1 1\n", "0 0 1\n", "0 4 1\n", "bad matrix\n"):
            with self.subTest(content=content):
                path = self.config()
                path.write_text(path.read_text().replace("FLAVORS=2", "FLAVORS=4")
                                .replace("[input]", '[input]\ninteraction_matrix="matrix.dat"'))
                (self.directory / "matrix.dat").write_text(content)
                self.validate_preserving_files(path)

    def test_cthyb_nonzero_retarded_kernel_two_iterations(self):
        path = self.config("hybridization", retarded=True)
        kernel = self.directory / "retarded.dat"
        kernel.write_text("".join(f"{tau:.17g} {0.01 * tau * (2 - tau):.17g} "
                                  f"{0.01 * (2 - 2 * tau):.17g}\n"
                                  for tau in (i / 8 for i in range(17))))
        original = path.read_bytes(), kernel.read_bytes()
        self.validate_preserving_files(path, valid=True)
        result = self.invoke(path)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual((path.read_bytes(), kernel.read_bytes()), original)
        self.assertEqual(list(self.work.iterdir()), [])
        self.assertEqual({item.name for item in self.directory.iterdir()},
                         {"run.toml", "retarded.dat", "results.h5"})
        self.verify("retarded")

    def test_retarded_validation_rejects_invalid_data_and_preserves_outputs(self):
        rows = [f"{i / 8:.17g} 0 0\n" for i in range(17)]
        valid = "".join(rows)
        cases = {"empty": "", "missing row": "".join(rows[:-1]),
                 "extra row": valid + rows[0], "extra text": valid + "extra\n",
                 "wrong tau grid": valid.replace("0.125", "0.2", 1),
                 "nonfinite value": valid.replace("0 0 0", "0 nan 0", 1),
                 "nonzero K at zero": valid.replace("0 0 0", "0 0.1 0", 1),
                 "negative K": valid.replace("0.125 0 0", "0.125 -0.1 0", 1)}
        for name, document in cases.items():
            with self.subTest(name=name):
                path = self.config("hybridization", retarded=True)
                (self.directory / "retarded.dat").write_text(document)
                (self.directory / "results.h5").write_bytes(b"existing scientific results")
                result = self.validate_preserving_files(path)
                self.assertIn("retarded", result.stderr.lower())

    def test_retarded_output_and_text_cannot_overwrite_the_kernel(self):
        path = self.config("hybridization", retarded=True)
        path.write_text(path.read_text().replace('results="results.h5"', 'results="retarded.dat"'))
        self.validate_preserving_files(path)
        path = self.config("hybridization", text=True, retarded=True)
        kernel = self.directory / "retarded.dat"
        alias = self.directory / "text" / "G_tau"
        try:
            alias.symlink_to(kernel)
        except OSError:
            self.skipTest("file symlinks unavailable")
        self.validate_preserving_files(path)

    def test_band_hopping_indices_match_spin_paired_band_count(self):
        path = self.config()
        document = path.read_text().replace("t=1.0", "t=1.0\nt0=1.0\nEPS_1=0.0\nEPSSQ_1=1.0")
        path.write_text(document)
        self.validate_preserving_files(path, valid=True)
        path.write_text(document.replace("t0=1.0", "t0=1.0\nt1=1.0"))
        result = self.validate_preserving_files(path)
        self.assertIn("t1", result.stderr)

    def test_custom_solver_selection_is_explicit(self):
        (self.directory / "custom.schema.toml").write_bytes(self.cthyb_schema.read_bytes())
        path = self.config("hybridization", loop="tau")
        builtin = path.read_text()
        schema = ("[input]", '[input]\nsolver_schema="custom.schema.toml"')
        delta = ("loop=", 'solver_input="delta"\nloop=')
        custom = builtin.replace('solver="hybridization"', f"solver={json.dumps(str(self.hybridization))}")
        cases = {"built-in with schema": builtin.replace(*schema),
                 "built-in with input kind": builtin.replace(*delta),
                 "custom without schema": custom.replace(*delta),
                 "custom without input kind": custom.replace(*schema)}
        for name, document in cases.items():
            with self.subTest(name=name):
                path.write_text(document)
                self.validate_preserving_files(path)
        path.write_text(custom.replace(*schema).replace(*delta))
        self.validate_preserving_files(path, valid=True)
        result = self.invoke(path)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.verify("sampled")

    def test_wide_external_sweep_counters_validate_without_sampling(self):
        path = self.config()
        path.write_text(path.read_text().replace("SWEEPS=32", "SWEEPS=1099511627776")
                        .replace("THERMALIZATION=4", "THERMALIZATION=1099511627776"))
        self.validate_preserving_files(path, valid=True)

    def test_unsupported_tau_solver_and_unsafe_dimensions_are_rejected(self):
        changes = [('loop="omega"', 'loop="tau"'),
                   ("N=16", "N=2147483646\n"),
                   ("NMATSUBARA=8", "NMATSUBARA=2147483647\n")]
        for old, new in changes:
            with self.subTest(new=new):
                path = self.config()
                path.write_text(path.read_text().replace(old, new).replace("FLAVORS=2", "FLAVORS=4"))
                self.validate_preserving_files(path)

    def test_frequency_measurements_cannot_overrun_the_grid(self):
        path = self.config()
        path.write_text(path.read_text().replace("NMATSUBARA=8", "NMATSUBARA=8\nNMATSUBARA_MEASUREMENTS=9"))
        result = self.validate_preserving_files(path)
        self.assertIn("NMATSUBARA_MEASUREMENTS", result.stderr)

    def test_text_output_cannot_follow_an_alias_to_the_run_file(self):
        path = self.config(text=True)
        alias = self.directory / "text" / "G_tau"
        try:
            alias.symlink_to(path)
        except OSError:
            self.skipTest("file symlinks unavailable")
        result = self.validate_preserving_files(path)
        self.assertIn("text output", result.stderr)

    def test_ctint_two_omega_iterations_and_analytic_reference(self):
        self.run_case("interaction", "exact")

    def test_text_output_is_explicit_and_relative(self):
        self.run_case("interaction", "exact", text=True)

    def test_invalid_configuration_preserves_existing_output(self):
        path = self.config()
        document = path.read_text(encoding="utf-8")
        output = self.directory / "results.h5"
        sentinel = b"existing scientific output"
        changes = [("BETA=2.0", "BETA=0.0"),
                   ("BETA=2.0", "BETA=2.0\nUNKNOWN_SCIENTIFIC_KEY=1"),
                   ("time_limit=0", "time_limit=-1"),
                   ("SYMMETRIZATION=true", "SYMMETRIZATION=true\nANTIFERROMAGNET=true"),
                   ('results="results.h5"', 'results="results.h5"\nfinal_tau="results.h5"')]
        for old, new in changes:
            with self.subTest(new=new):
                path.write_text(document.replace(old, new), encoding="utf-8")
                output.write_bytes(sentinel)
                result = self.invoke("--validate", path)
                self.assertNotEqual(result.returncode, 0)
                self.assertEqual(output.read_bytes(), sentinel)
                self.assertEqual(list(self.work.iterdir()), [])

    def test_output_cannot_replace_run_file(self):
        path = self.config()
        path.write_text(path.read_text(encoding="utf-8").replace('results="results.h5"', 'results="run.toml"'),
                        encoding="utf-8")
        original = path.read_bytes()
        result = self.invoke("--validate", path)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(path.read_bytes(), original)

    def test_legacy_run_and_invalid_cli_rejected(self):
        path = self.directory / "legacy.parm"
        path.write_text("BETA=2; U=0; N=16; SOLVER=hybridization;", encoding="utf-8")
        for arguments in [(), (path,), ("--unknown",), ("missing.toml",), ("a.toml", "b.toml")]:
            with self.subTest(arguments=arguments):
                self.assertNotEqual(self.invoke(*arguments).returncode, 0)
        self.assertEqual(self.invoke("--help", "missing.toml").returncode, 0)
        self.assertEqual({p.name for p in self.directory.iterdir()}, {"legacy.parm"})
        self.assertEqual(list(self.work.iterdir()), [])

    def validate_preserving_files(self, path, valid=False):
        before = {item.relative_to(self.directory): item.read_bytes()
                  for item in self.directory.rglob("*") if item.is_file()}
        result = self.invoke("--validate", path)
        self.assertEqual({item.relative_to(self.directory): item.read_bytes()
                          for item in self.directory.rglob("*") if item.is_file()}, before)
        self.assertEqual(list(self.work.iterdir()), [])
        if valid:
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        else:
            self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
        return result

    def test_dos_validation_rejects_unsafe_grids_and_weights(self):
        cases = {
            "empty": "",
            "one point": "0 1\n",
            "even point count": "-1 1\n1 1\n",
            "truncated pair": "-1 1\n0 1\n1\n",
            "trailing garbage": "-1 1\n0 1\n1 1\nextra\n",
            "nonfinite energy": "-1 1\n0 1\ninf 1\n",
            "nonfinite density": "-1 1\n0 nan\n1 1\n",
            "negative density": "-1 1\n0 -1\n1 1\n",
            "zero weight": "-1 0\n0 0\n1 0\n",
            "descending": "1 1\n0 1\n-1 1\n",
            "duplicate energy": "-1 1\n-1 1\n1 1\n",
            "unequal spacing": "-1 1\n0.25 1\n1 1\n",
        }
        for name, document in cases.items():
            with self.subTest(name=name):
                path = self.config()
                path.write_text(path.read_text().replace("[input]", '[input]\ndos="dos.dat"'))
                (self.directory / "dos.dat").write_text(document)
                (self.directory / "results.h5").write_bytes(b"existing scientific results")
                result = self.validate_preserving_files(path)
                self.assertIn("dos", result.stderr.lower())

    def test_default_general_transform_accepts_a_valid_dos(self):
        path = self.config()
        path.write_text(path.read_text().replace("[input]", '[input]\ndos="dos.dat"'))
        # Zero endpoints are valid; the integrated weight must be positive.
        (self.directory / "dos.dat").write_text("-1 0\n0 1\n1 0\n")
        (self.directory / "results.h5").write_bytes(b"existing scientific results")
        self.validate_preserving_files(path, valid=True)

    def test_initial_green_text_validation(self):
        for section, loop, solver in [("initial_omega", "omega", "interaction"),
                                      ("initial_tau", "tau", "hybridization")]:
            with self.subTest(section=section):
                rows = ([f"{i} (0,{-2.0 / ((2*i+1)*math.pi):.17g}) "
                         f"(0,{-2.0 / ((2*i+1)*math.pi):.17g})\n" for i in range(8)]
                        if section == "initial_omega" else
                        [f"{i} -0.5 -0.5\n" for i in range(17)])
                valid = "".join(rows)
                path = self.config(solver, loop)
                path.write_text(path.read_text().replace("[input]",
                                                        f'[input]\n{section}="initial.dat"'))
                (self.directory / "initial.dat").write_text(valid)
                (self.directory / "results.h5").write_bytes(b"existing scientific results")
                self.validate_preserving_files(path, valid=True)
                cases = ["", "".join(rows[:-1]), valid + rows[0], valid + "extra\n",
                         valid.replace("0 ", "nan ", 1),
                         valid.replace("(0,", "(nan," , 1) if section == "initial_omega" else
                         valid.replace("-0.5", "nan", 1)]
                for document in cases:
                    with self.subTest(document=document[:40]):
                        (self.directory / "initial.dat").write_text(document)
                        result = self.validate_preserving_files(path)
                        self.assertIn("green", result.stderr.lower())

    def test_square_lattice_retains_tprime_and_grid_settings(self):
        path = self.config()
        document = path.read_text().replace('t=1.0', 't=1.0\nTWODBS="square"\ntprime=0.25\nL=8')
        path.write_text(document)
        (self.directory / "results.h5").write_bytes(b"existing scientific results")
        self.validate_preserving_files(path, valid=True)
        for old, new in [("L=8", "L=0"), ("L=8", "L=8.5"),
                         ("tprime=0.25", "tprime=nan")]:
            with self.subTest(new=new):
                path.write_text(document.replace(old, new))
                self.validate_preserving_files(path)

    def test_child_failure_propagates_without_replacing_results(self):
        path = self.config()
        # The native helper accepts this TOML protocol and deliberately exits 7.
        (self.directory / "helper.schema.toml").write_text('''application="dmft-contract"
schema_version=1
[parameters.N]
type="int64"
required=true
[parameters.NMATSUBARA]
type="int64"
required=true
[parameters.FLAVORS]
type="int64"
required=true
[parameters.MU]
type="float64"
required=true
[parameters.MODE]
type="string"
default="fail"
[input.g0]
type="path"
required=true
[output.results]
type="path"
required=true
[execution.seed]
type="int64"
required=true
''', encoding="utf-8")
        text = path.read_text(encoding="utf-8").replace(
            'solver="interaction"', f'solver={json.dumps(str(self.helper))}\nsolver_input="g0"')
        text = text.replace("ALPHA=-0.01\nMEASUREMENT_PERIOD=1\n", "")
        text = text.replace("SWEEPS=32\nTHERMALIZATION=4\n", "")
        text = text.replace("[input]", '[input]\nsolver_schema="helper.schema.toml"')
        path.write_text(text, encoding="utf-8")
        output = self.directory / "results.h5"
        output.write_bytes(b"existing scientific output")
        result = self.invoke(path)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("7", result.stderr)
        self.assertIn("solver", result.stderr.lower())
        self.assertEqual(output.read_bytes(), b"existing scientific output")
        self.assertEqual(list(self.work.iterdir()), [])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("dmft", type=Path)
    parser.add_argument("hybridization", type=Path)
    parser.add_argument("helper", type=Path)
    arguments = parser.parse_args()
    for name, path in vars(arguments).items():
        setattr(DMFTCLIContract, name, path.resolve(strict=True))
    DMFTCLIContract.cthyb_schema = Path(__file__).resolve().parents[1] / "hybridization/schema/cthyb.toml"
    result = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(DMFTCLIContract))
    raise SystemExit(0 if result.wasSuccessful() else 1)
