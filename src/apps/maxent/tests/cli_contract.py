"""Exercise TOML input and fail-before-output behavior with the real binary."""
import argparse
from pathlib import Path
import subprocess
import tempfile
import unittest

class MaxEntCLIContract(unittest.TestCase):
    def setUp(self):
        temporary=tempfile.TemporaryDirectory(prefix="alps-maxent-cli-")
        self.addCleanup(temporary.cleanup)
        self.directory=Path(temporary.name)
    def run_cli(self,*args):
        before={p.name:p.read_bytes() for p in self.directory.iterdir() if p.is_file()}
        result=subprocess.run([str(self.executable),*map(str,args)],cwd=self.directory,
                              capture_output=True,text=True,timeout=30)
        self.assertEqual({p.name:p.read_bytes() for p in self.directory.iterdir() if p.is_file()},before)
        return result
    def test_help_without_loading_input(self):
        result=self.run_cli('--help','missing.toml')
        self.assertEqual(result.returncode,0)
        self.assertIn('--validate',result.stdout)
        self.assertEqual(result.stderr,'')
    def test_rejected_cli(self):
        for args in [(),('--bad-option',),('missing.toml',)]:
            with self.subTest(args=args):
                result=self.run_cli(*args)
                self.assertEqual(result.returncode,1)
                self.assertIn('maxent:',result.stderr)
    def config(self,extra=''):
        path=self.directory/'run.toml'
        path.write_text('''format_version=1
application="maxent"
schema_version=1
[parameters]
BETA=2.0
NFREQ=20
OMEGA_MAX=4.0
'''+extra+'''
[input]
values=[-0.5,-0.3,-0.3,-0.5]
errors=[0.01,0.01,0.01,0.01]
[output]
results="result.h5"
''')
        return path
    def test_validate_has_no_output(self):
        result=self.run_cli('--validate',self.config())
        self.assertEqual(result.returncode,0,result.stderr)
    def test_unknown_key_does_not_touch_existing_output(self):
        (self.directory/'result.h5').write_bytes(b'preserve existing result')
        result=self.run_cli(self.config('NFRQ=20'))
        self.assertEqual(result.returncode,1)
        self.assertIn('parameters.NFRQ',result.stderr)
    def test_invalid_scientific_settings_do_not_touch_outputs(self):
        for setting in ['T=1.0', 'NFREQ=3', 'DEFAULT_MODEL="gaussian"',
                        'FREQUENCY_GRID="log"\nNFREQ=21', 'KERNEL="anomalous"']:
            with self.subTest(setting=setting):
                path=self.config()
                text=path.read_text()
                # Replace existing fields, append optional fields in the parameter section.
                for line in setting.splitlines():
                    key=line.split('=')[0]
                    import re
                    if re.search(r'^'+key+r'=.*$',text,re.M):
                        text=re.sub(r'^'+key+r'=.*$',lambda _:line,text,flags=re.M)
                    else:
                        text=text.replace('[input]',line+'\n[input]')
                path.write_text(text)
                (self.directory/'result.h5').write_bytes(b'preserve existing result')
                result=self.run_cli('--validate',path)
                self.assertEqual(result.returncode,1,result.stdout)

    def test_legacy_input_rejected(self):
        path=self.directory/'old.parm'
        path.write_text('BETA=2; NFREQ=20;')
        self.assertEqual(self.run_cli(path).returncode,1)

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('maxent',type=Path)
    args=parser.parse_args()
    MaxEntCLIContract.executable=args.maxent.resolve(strict=True)
    result=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(MaxEntCLIContract))
    raise SystemExit(0 if result.wasSuccessful() else 1)
