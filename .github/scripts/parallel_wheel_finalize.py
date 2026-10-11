"""Prototype: find the separately installed MPICH wheel in a Unix venv.

MPICH installs libmpi into the environment's lib directory. Do not bundle a
second copy into pyalps: its launcher and mpi4py must use the same library.
This intentionally tests the standard venv layout before proposing a release.
"""
import os
from pathlib import Path
import subprocess
import sys
import tempfile

for wheel in Path(sys.argv[1]).glob('*.whl'):
    with tempfile.TemporaryDirectory() as temporary:
        subprocess.run([sys.executable, '-m', 'wheel', 'unpack', str(wheel), '-d', temporary], check=True)
        tree, = Path(temporary).iterdir()
        # tree corresponds to <prefix>/lib/pythonX.Y/site-packages.
        prefix_lib = tree / '../..'
        for binary in tree.rglob('*'):
            if not binary.is_file():
                continue
            if sys.platform == 'darwin':
                result = subprocess.run(['otool', '-L', str(binary)], capture_output=True, text=True)
                changes = []
                for line in result.stdout.splitlines()[1:]:
                    dependency = line.strip().split(' (', 1)[0]
                    if Path(dependency).name == 'libmpi.12.dylib':
                        relative = os.path.relpath(prefix_lib / 'libmpi.12.dylib', binary.parent)
                        changes += ['-change', dependency, '@loader_path/' + relative]
                if changes:
                    subprocess.run(['install_name_tool', *changes, str(binary)], check=True)
                    subprocess.run(['codesign', '--force', '--sign', '-', str(binary)], check=True)
                continue
            with binary.open('rb') as stream:
                if stream.read(4) != b'\x7fELF':
                    continue
            needed = subprocess.check_output(['patchelf', '--print-needed', str(binary)], text=True)
            if 'libmpi.so.12' in needed.splitlines():
                relative = os.path.relpath(prefix_lib, binary.parent)
                subprocess.run(['patchelf', '--add-rpath', '$ORIGIN/' + relative, str(binary)], check=True)
        subprocess.run([sys.executable, '-m', 'wheel', 'pack', str(tree), '-d', str(wheel.parent)], check=True)
