"""Prepare pinned CI dependencies using official sources and reusable caches.

CI builds selected Boost releases when the OS package is too old. HDF5 uses
prebuilt conda-forge packages on glibc and a minimal C-only build on musl.
Windows uses vcpkg's binary cache and builds any missing packages.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import subprocess
import tarfile
import tempfile
import urllib.request

ROOT = Path(__file__).resolve().parents[2]

HEADER_ONLY = {
    'tomlplusplus': ('marzer/tomlplusplus', [
        '-DTOMLPLUSPLUS_BUILD_TESTS=OFF', '-DTOMLPLUSPLUS_BUILD_EXAMPLES=OFF',
        '-DTOMLPLUSPLUS_BUILD_MODULES=OFF',
    ]),
    'highfive': ('highfive-devs/highfive', [
        '-DHIGHFIVE_FIND_HDF5=OFF', '-DHIGHFIVE_UNIT_TESTS=OFF',
        '-DHIGHFIVE_EXAMPLES=OFF', '-DHIGHFIVE_BUILD_DOCS=OFF',
    ]),
}


def download_verified(url: str, archive: Path, checksum: str):
    digest = hashlib.sha256()
    with urllib.request.urlopen(url, timeout=120) as response, archive.open('wb') as stream:
        while chunk := response.read(1024 * 1024):
            digest.update(chunk)
            stream.write(chunk)
    if digest.hexdigest() != checksum:
        raise ValueError(f'Checksum mismatch for {archive.name}')


def dependency_cache(destination: Path, name: str, fingerprint: str):
    stamp = destination / f'.alps-{name}-sha256'
    if stamp.is_file() and stamp.read_text() == fingerprint:
        return stamp, True
    if destination.exists() and any(destination.iterdir()):
        raise ValueError(f'Refusing to replace an existing dependency: {destination}')
    destination.parent.mkdir(parents=True, exist_ok=True)
    return stamp, False


def prepare_boost(version: str, destination: Path):
    checksums = json.loads((ROOT / '.github/ci-matrix.json').read_text())['boost']
    if version not in checksums:
        raise ValueError(f'Unpinned Boost version: {version}')
    fingerprint = checksums[version] + os.environ.get('CXX', 'c++') + os.environ.get('ALPS_BOOST_MPI', 'ON')
    stamp, cached = dependency_cache(destination, 'boost', fingerprint)
    if cached:
        return
    with tempfile.TemporaryDirectory(prefix='.alps-boost-', dir=destination.parent) as temporary:
        archive = Path(temporary) / 'boost.tar.gz'
        url = f'https://archives.boost.io/release/{version}/source/boost_{version.replace(".", "_")}.tar.gz'
        download_verified(url, archive, checksums[version])
        source = Path(temporary) / 'source'
        with tarfile.open(archive) as package:
            package.extractall(source, filter='data')
        extracted, = source.iterdir()
        subprocess.run(['bash', str(ROOT / '.github/scripts/build-boost.sh'), str(extracted), str(destination.resolve())], check=True)
    stamp.write_text(fingerprint)


def prepare_header_only(name: str, destination: Path):
    """Install a pinned header-only CMake package without compiling dependencies."""
    repository, options = HEADER_ONLY[name]
    pin = json.loads((ROOT / '.github/dependencies.json').read_text())[name]
    stamp, cached = dependency_cache(destination, name, pin['sha256'])
    if cached:
        return
    with tempfile.TemporaryDirectory(prefix=f'.alps-{name}-', dir=destination.parent) as temporary:
        archive = Path(temporary) / f'{name}.tar.gz'
        url = f'https://codeload.github.com/{repository}/tar.gz/refs/tags/v{pin["version"]}'
        download_verified(url, archive, pin['sha256'])
        source = Path(temporary) / 'source'
        with tarfile.open(archive) as package:
            package.extractall(source, filter='data')
        extracted, = source.iterdir()
        build = Path(temporary) / 'build'
        subprocess.run(['cmake', '-S', str(extracted), '-B', str(build),
                        *options,
                        f'-DCMAKE_INSTALL_PREFIX={destination}'], check=True)
        subprocess.run(['cmake', '--install', str(build)], check=True)
    stamp.write_text(pin['sha256'])


def prepare_hdf5(destination: Path, *, musl=False):
    """Install one HDF5 C provider without adding compilers or a Python environment."""
    if platform.system() != 'Linux':
        raise ValueError('Use Homebrew or vcpkg for HDF5 on non-Linux platforms')
    pins = json.loads((ROOT / '.github/dependencies.json').read_text())
    pin = pins['hdf5']
    if musl:
        fingerprint = pin['sha256'] + platform.machine() + os.environ.get('CC', 'cc')
    else:
        architecture = {'x86_64': 'linux-64', 'aarch64': 'linux-aarch64'}[platform.machine()]
        fingerprint = pin['version'] + pin['conda_builds'][architecture] + pins['micromamba']['release']
    stamp, cached = dependency_cache(destination, 'hdf5', fingerprint)
    if cached:
        return
    with tempfile.TemporaryDirectory(prefix='.alps-hdf5-', dir=destination.parent) as temporary:
        temporary = Path(temporary)
        if musl:
            archive = temporary / 'hdf5.tar.gz'
            url = f'https://github.com/HDFGroup/hdf5/releases/download/{pin["version"]}/hdf5-{pin["version"]}.tar.gz'
            download_verified(url, archive, pin['sha256'])
            with tarfile.open(archive) as package:
                package.extractall(temporary / 'source', filter='data')
            source, = (temporary / 'source').iterdir()
            build = temporary / 'build'
            subprocess.run(['cmake', '-S', str(source), '-B', str(build), '-G', 'Ninja',
                            '-DCMAKE_BUILD_TYPE=Release', '-DBUILD_SHARED_LIBS=ON',
                            '-DBUILD_STATIC_LIBS=OFF', '-DBUILD_TESTING=OFF',
                            '-DHDF5_BUILD_TOOLS=OFF', '-DHDF5_BUILD_EXAMPLES=OFF',
                            '-DHDF5_BUILD_CPP_LIB=OFF', '-DHDF5_BUILD_FORTRAN=OFF',
                            '-DHDF5_BUILD_JAVA=OFF', '-DHDF5_BUILD_HL_LIB=OFF',
                            '-DHDF5_ENABLE_ZLIB_SUPPORT=ON', '-DHDF5_ENABLE_SZIP_SUPPORT=OFF',
                            f'-DCMAKE_INSTALL_PREFIX={destination}'], check=True)
            subprocess.run(['cmake', '--build', str(build), '--parallel'], check=True)
            subprocess.run(['cmake', '--install', str(build)], check=True)
        else:
            tool_pin = pins['micromamba']
            archive = temporary / 'micromamba.tar.bz2'
            url = f'https://github.com/mamba-org/micromamba-releases/releases/download/{tool_pin["release"]}/micromamba-{architecture}.tar.bz2'
            download_verified(url, archive, tool_pin['sha256'][architecture])
            tool = temporary / 'micromamba'
            with tarfile.open(archive) as package:
                tool.write_bytes(package.extractfile('bin/micromamba').read())
            tool.chmod(0o755)
            subprocess.run([str(tool), '--no-rc', '--root-prefix', str(temporary / 'cache'),
                            'create', '--yes', '--prefix', str(destination),
                            '--override-channels', '--channel', 'conda-forge',
                            '--strict-channel-priority',
                            f'hdf5={pin["version"]}={pin["conda_builds"][architecture]}'],
                           check=True, env={**os.environ, 'CONDA_PKGS_DIRS': str(temporary / 'packages')})
    stamp.write_text(fingerprint)


def prepare_windows(architecture: str, destination: Path, *, environment=None):
    if architecture not in {'x64', 'arm64'}:
        raise ValueError(f'Unsupported architecture: {architecture}')
    baseline = json.loads((ROOT / 'cmake/vcpkg/vcpkg.json').read_text())['builtin-baseline']
    if not re.fullmatch('[a-f0-9]{40}', baseline):
        raise ValueError('Invalid vcpkg baseline')
    tool = destination / 'vcpkg.exe'
    if not (destination / '.git').exists():
        if destination.exists() and any(destination.iterdir()):
            raise ValueError(f'Refusing to replace an existing directory: {destination}')
        subprocess.run(['git', 'clone', '--filter=blob:none', '--no-checkout',
                        'https://github.com/microsoft/vcpkg.git', str(destination)], check=True, env=environment)
    subprocess.run(['git', '-C', str(destination), 'checkout', '--detach', baseline], check=True, env=environment)
    subprocess.run([str(destination / 'bootstrap-vcpkg.bat'), '-disableMetrics'], check=True, env=environment)
    subprocess.run([str(tool), 'install', f'--triplet={architecture}-windows',
                    f'--x-manifest-root={ROOT / "cmake/vcpkg"}', f'--x-install-root={destination / "installed"}',
                    '--x-feature=tests',
                    *([f'--overlay-ports={ROOT / "cmake/vcpkg-arm64-overlay"}'] if architecture == 'arm64' else [])],
                   check=True, env=environment)
    return destination


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('package')
    parser.add_argument('destination', type=Path)
    args = parser.parse_args()
    if args.package in HEADER_ONLY:
        prepare_header_only(args.package, args.destination.resolve())
    elif args.package in {'hdf5', 'hdf5-musl'}:
        prepare_hdf5(args.destination.resolve(), musl=args.package == 'hdf5-musl')
    elif args.package.startswith('windows-'):
        prepare_windows(args.package.removeprefix('windows-'), args.destination.resolve())
    else:
        match = re.fullmatch(r'boost-(\d+\.\d+\.\d+)(?:-[a-z0-9-]+)?', args.package)
        if not match:
            parser.error('Expected boost-VERSION, tomlplusplus, highfive, hdf5, hdf5-musl or windows-{x64,arm64}')
        prepare_boost(match[1], args.destination.resolve())


if __name__ == '__main__':
    main()
