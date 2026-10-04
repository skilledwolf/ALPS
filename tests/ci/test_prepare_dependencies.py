"""Check source integrity and architecture selection without downloading packages."""
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import tarfile

import pytest

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location('prepare_dependencies', ROOT / '.github/scripts/prepare_dependencies.py')
dependencies = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(dependencies)


def test_boost_checksum_failure_never_extracts_or_builds(tmp_path, monkeypatch):
    monkeypatch.setattr(dependencies.urllib.request, 'urlopen', lambda *a, **kw: io.BytesIO(b'corrupt archive'))
    def forbidden(*a, **kw):
        pytest.fail('An unverified dependency must not be extracted or built')
    monkeypatch.setattr(dependencies.tarfile, 'open', forbidden)
    monkeypatch.setattr(dependencies.subprocess, 'run', forbidden)
    with pytest.raises(ValueError, match='Checksum mismatch'):
        dependencies.prepare_boost('1.91.0', tmp_path / 'installed')
    assert not (tmp_path / 'installed').exists()


def test_matching_boost_cache_avoids_network(tmp_path, monkeypatch):
    monkeypatch.setenv('CXX', 'clang++')
    monkeypatch.setenv('ALPS_BOOST_MPI', 'OFF')
    checksum = json.loads((ROOT / '.github/ci-matrix.json').read_text())['boost']['1.91.0']
    (tmp_path / '.alps-boost-sha256').write_text(checksum + 'clang++OFF')
    monkeypatch.setattr(dependencies.urllib.request, 'urlopen', lambda *a, **kw: pytest.fail('Cached dependency downloaded again'))
    dependencies.prepare_boost('1.91.0', tmp_path)


@pytest.mark.parametrize('architecture', ['x64', 'arm64'])
def test_windows_uses_pinned_baseline_and_matching_triplet(tmp_path, monkeypatch, architecture):
    commands = []
    monkeypatch.setattr(dependencies.subprocess, 'run', lambda args, **kw: commands.append(args))
    dependencies.prepare_windows(architecture, tmp_path)
    baseline = json.loads((ROOT / 'cmake/vcpkg/vcpkg.json').read_text())['builtin-baseline']
    assert ['git', '-C', str(tmp_path), 'checkout', '--detach', baseline] in commands
    install = commands[-1]
    assert f'--triplet={architecture}-windows' in install
    assert f'--x-manifest-root={ROOT / "cmake/vcpkg"}' in install
    assert '--x-feature=tests' in install
    assert any(arg.startswith('--overlay-ports=') for arg in install) == (architecture == 'arm64')


@pytest.mark.parametrize('name', ['tomlplusplus', 'highfive'])
def test_header_only_checksum_failure_never_installs(tmp_path, monkeypatch, name):
    monkeypatch.setattr(dependencies.urllib.request, 'urlopen', lambda *a, **kw: io.BytesIO(b'corrupt'))
    monkeypatch.setattr(dependencies.tarfile, 'open', lambda *a, **kw: pytest.fail('Unverified package extracted'))
    monkeypatch.setattr(dependencies.subprocess, 'run', lambda *a, **kw: pytest.fail('Unverified package installed'))
    with pytest.raises(ValueError, match='Checksum mismatch'):
        dependencies.prepare_header_only(name, tmp_path / 'installed')
    assert not (tmp_path / 'installed').exists()


@pytest.mark.parametrize('name', ['tomlplusplus', 'highfive'])
def test_matching_header_only_cache_avoids_network(tmp_path, monkeypatch, name):
    pin = json.loads((ROOT / '.github/dependencies.json').read_text())[name]
    (tmp_path / f'.alps-{name}-sha256').write_text(pin['sha256'])
    monkeypatch.setattr(dependencies.urllib.request, 'urlopen', lambda *a, **kw: pytest.fail('Cached package downloaded again'))
    dependencies.prepare_header_only(name, tmp_path)


@pytest.mark.parametrize('name', ['tomlplusplus', 'highfive'])
def test_header_only_provider_installs_without_building(tmp_path, monkeypatch, name):
    payload = io.BytesIO()
    with tarfile.open(fileobj=payload, mode='w:gz') as archive:
        entry = tarfile.TarInfo('source/CMakeLists.txt')
        contents = b'project(header_only)'
        entry.size = len(contents)
        archive.addfile(entry, io.BytesIO(contents))
    checksum = hashlib.sha256(payload.getvalue()).hexdigest()
    pins = tmp_path / 'repository/.github'
    pins.mkdir(parents=True)
    (pins / 'dependencies.json').write_text(json.dumps({name: {'version': '3.3.0', 'sha256': checksum}}))
    monkeypatch.setattr(dependencies, 'ROOT', pins.parent)
    monkeypatch.setattr(dependencies.urllib.request, 'urlopen', lambda *a, **kw: io.BytesIO(payload.getvalue()))
    installed = tmp_path / 'installed'
    commands = []
    def run(args, **kw):
        commands.append(args)
        if '--install' in args:
            installed.mkdir()
    monkeypatch.setattr(dependencies.subprocess, 'run', run)
    dependencies.prepare_header_only(name, installed)
    assert len(commands) == 2
    assert commands[0][:2] == ['cmake', '-S']
    assert all(option in commands[0] for option in dependencies.HEADER_ONLY[name][1])
    assert commands[1][:2] == ['cmake', '--install']
    assert (installed / f'.alps-{name}-sha256').read_text() == checksum


@pytest.mark.parametrize('musl', [False, True])
def test_hdf5_checksum_failure_never_extracts_or_installs(tmp_path, monkeypatch, musl):
    monkeypatch.setattr(dependencies.platform, 'system', lambda: 'Linux')
    monkeypatch.setattr(dependencies.platform, 'machine', lambda: 'x86_64')
    monkeypatch.setattr(dependencies.urllib.request, 'urlopen', lambda *a, **kw: io.BytesIO(b'corrupt'))
    monkeypatch.setattr(dependencies.tarfile, 'open', lambda *a, **kw: pytest.fail('Unverified archive extracted'))
    monkeypatch.setattr(dependencies.subprocess, 'run', lambda *a, **kw: pytest.fail('Unverified provider installed'))
    with pytest.raises(ValueError, match='Checksum mismatch'):
        dependencies.prepare_hdf5(tmp_path / 'installed', musl=musl)
    assert not (tmp_path / 'installed').exists()


@pytest.mark.parametrize('machine,architecture', [('x86_64', 'linux-64'), ('aarch64', 'linux-aarch64')])
def test_glibc_hdf5_provider_uses_pinned_serial_build_and_scoped_cache(tmp_path, monkeypatch, machine, architecture):
    monkeypatch.setattr(dependencies.platform, 'system', lambda: 'Linux')
    monkeypatch.setattr(dependencies.platform, 'machine', lambda: machine)
    payload = io.BytesIO()
    with tarfile.open(fileobj=payload, mode='w:bz2') as archive:
        entry = tarfile.TarInfo('bin/micromamba')
        contents = b'#!/bin/sh\n'
        entry.size = len(contents)
        archive.addfile(entry, io.BytesIO(contents))
    pins = json.loads((ROOT / '.github/dependencies.json').read_text())
    pins['micromamba']['sha256'][architecture] = hashlib.sha256(payload.getvalue()).hexdigest()
    repository = tmp_path / 'repository/.github'
    repository.mkdir(parents=True)
    (repository / 'dependencies.json').write_text(json.dumps(pins))
    monkeypatch.setattr(dependencies, 'ROOT', repository.parent)
    urls = []
    def download(url, **kw):
        urls.append(url)
        return io.BytesIO(payload.getvalue())
    monkeypatch.setattr(dependencies.urllib.request, 'urlopen', download)
    installed = tmp_path / 'installed'
    commands = []
    def run(args, **kw):
        commands.append((args, kw))
        installed.mkdir()
    monkeypatch.setattr(dependencies.subprocess, 'run', run)
    dependencies.prepare_hdf5(installed)
    assert len(commands) == 1
    command, options = commands[0]
    assert '--no-rc' in command and '--override-channels' in command
    assert command[-1] == f'hdf5=2.2.0={pins["hdf5"]["conda_builds"][architecture]}'
    assert '--prefix' in command and str(installed) in command
    assert '--channel' in command and 'conda-forge' in command
    assert any(f'micromamba-{architecture}.tar.bz2' in url for url in urls)
    assert not Path(options['env']['CONDA_PKGS_DIRS']).exists()
    monkeypatch.setattr(dependencies.urllib.request, 'urlopen', lambda *a, **kw: pytest.fail('Cached provider downloaded again'))
    dependencies.prepare_hdf5(installed)
    assert len(commands) == 1


def test_musl_hdf5_provider_only_builds_shared_c_library(tmp_path, monkeypatch):
    monkeypatch.setattr(dependencies.platform, 'system', lambda: 'Linux')
    monkeypatch.setattr(dependencies.platform, 'machine', lambda: 'x86_64')
    payload = io.BytesIO()
    with tarfile.open(fileobj=payload, mode='w:gz') as archive:
        entry = tarfile.TarInfo('source/CMakeLists.txt')
        contents = b'project(hdf5)'
        entry.size = len(contents)
        archive.addfile(entry, io.BytesIO(contents))
    pins = json.loads((ROOT / '.github/dependencies.json').read_text())
    pins['hdf5']['sha256'] = hashlib.sha256(payload.getvalue()).hexdigest()
    repository = tmp_path / 'repository/.github'
    repository.mkdir(parents=True)
    (repository / 'dependencies.json').write_text(json.dumps(pins))
    monkeypatch.setattr(dependencies, 'ROOT', repository.parent)
    monkeypatch.setattr(dependencies.urllib.request, 'urlopen', lambda *a, **kw: io.BytesIO(payload.getvalue()))
    installed = tmp_path / 'installed'
    commands = []
    def run(args, **kw):
        commands.append(args)
        if '--install' in args:
            installed.mkdir()
    monkeypatch.setattr(dependencies.subprocess, 'run', run)
    dependencies.prepare_hdf5(installed, musl=True)
    assert len(commands) == 3
    configure = commands[0]
    assert '-DBUILD_SHARED_LIBS=ON' in configure
    for option in ['BUILD_STATIC_LIBS', 'BUILD_TESTING', 'HDF5_BUILD_TOOLS', 'HDF5_BUILD_EXAMPLES',
                   'HDF5_BUILD_CPP_LIB', 'HDF5_BUILD_FORTRAN', 'HDF5_BUILD_JAVA', 'HDF5_BUILD_HL_LIB']:
        assert f'-D{option}=OFF' in configure
    assert '-DHDF5_ENABLE_ZLIB_SUPPORT=ON' in configure
    assert '-DHDF5_ENABLE_SZIP_SUPPORT=OFF' in configure
    assert commands[1][1] == '--build' and commands[2][1] == '--install'
    monkeypatch.setattr(dependencies.urllib.request, 'urlopen', lambda *a, **kw: pytest.fail('Cached provider downloaded again'))
    dependencies.prepare_hdf5(installed, musl=True)
    assert len(commands) == 3


@pytest.mark.parametrize('name', ['boost', 'highfive', 'hdf5'])
def test_dependency_cache_preserves_unrelated_files(tmp_path, name):
    user_file = tmp_path / 'notes.txt'
    user_file.write_text('keep this')
    with pytest.raises(ValueError, match='Refusing to replace'):
        dependencies.dependency_cache(tmp_path, name, 'different')
    assert user_file.read_text() == 'keep this'
