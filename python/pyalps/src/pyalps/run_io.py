# Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
"""Write and execute TOML runs and explicit multi-run job manifests.

The application supplies its schema separately. Arrays are parameter values;
each entry passed to ``write_run_files`` is one run, never an inferred sweep.
"""

from collections.abc import Mapping
import json
import operator
import os
from pathlib import Path
import subprocess
import tempfile
import tomllib

from . import run_config

_SECTIONS = ("parameters", "input", "output", "execution")


def _quoted(value):
    # TOML shares JSON's basic-string escapes but forbids surrogate escapes
    # and unescaped DEL. Keep complete Unicode characters in the UTF-8 file.
    return json.dumps(value, ensure_ascii=False).replace("\x7f", "\\u007f")


def _section(values):
    if values is None:
        return {}
    if not isinstance(values, Mapping):
        raise TypeError("run sections must be mappings")
    if any(not isinstance(key, str) for key in values):
        raise TypeError("run keys must be strings")
    return {key: os.fspath(value) if isinstance(value, os.PathLike) else value
            for key, value in values.items()}


def _document(schema, sections):
    supplied = {name: _section(sections.get(name)) for name in _SECTIONS}
    explicit = run_config.RunConfiguration()
    if schema is None:
        # The application validates the run before execution.
        for name in _SECTIONS:
            target = getattr(explicit, name)
            for key, value in supplied[name].items():
                target[key] = value
        return run_config.format(explicit), explicit
    # Preserve relative paths in the document. The application's native loader
    # resolves them against this file, rather than the writer's working directory.
    resolved = run_config.resolve(schema, **supplied, base_directory="")
    for name in _SECTIONS:
        source, target = getattr(resolved, name), getattr(explicit, name)
        for key in supplied[name]:
            target[key] = source[key]
    return run_config.format(explicit), resolved


def _check_output_paths(runs, targets):
    configurations = {path.resolve() for path in targets}
    owners = {}
    for index, (path, resolved) in enumerate(runs):
        for key in ("results", "checkpoint"):
            if key not in resolved.output:
                continue
            output = Path(resolved.output[key])
            if not output.is_absolute():
                output = path.parent / output
            output = output.resolve()
            if output in configurations:
                raise ValueError(f"output.{key} would overwrite a run or job file: {output}")
            if output in owners and owners[output] != index:
                raise ValueError(f"multiple runs would write output file: {output}")
            owners[output] = index


def _publish(documents, overwrite):
    """Write complete files, refusing existing destinations unless overwriting."""
    targets = [path for path, _ in documents]
    if len(set(targets)) != len(targets):
        raise ValueError("run and manifest filenames must be distinct")
    if not overwrite:
        for path in targets:
            if path.exists() or path.is_symlink():
                raise FileExistsError(path)
    created = []
    try:
        for path, document in documents:
            path.parent.mkdir(parents=True, exist_ok=True)
            if overwrite:
                descriptor, temporary = tempfile.mkstemp(dir=path.parent, prefix=".alps-run-")
                try:
                    with os.fdopen(descriptor, "w", encoding="utf-8", newline="\n") as stream:
                        stream.write(document)
                    os.replace(temporary, path)
                finally:
                    Path(temporary).unlink(missing_ok=True)
            else:
                with open(path, "x", encoding="utf-8", newline="\n") as stream:
                    created.append(path)
                    stream.write(document)
    except BaseException:
        for path in created:
            path.unlink(missing_ok=True)
        raise


def write_run_file(filename, schema=None, *, parameters=None, input=None, output=None,
                   execution=None, overwrite=False):
    """Write one run and return its path, omitting schema-supplied defaults.

    ``schema`` is the application-owned TOML schema text. When supplied, the
    native validator used by C++ applications validates all four run sections;
    otherwise the application validates the run when it is executed.
    """
    path = Path(filename)
    document, resolved = _document(schema, dict(parameters=parameters, input=input,
                                                output=output, execution=execution))
    _check_output_paths([(path, resolved)], [path])
    _publish([(path, document)], overwrite)
    return path


def write_run_files(prefix, runs, schema=None, *, baseseed=None, overwrite=False):
    """Write explicit runs and return their ``<prefix>.job.toml`` manifest.

    Each run maps section names to mappings. If ``baseseed`` is supplied,
    missing ``execution.seed`` values receive ``baseseed + run_index``;
    explicit seeds and the caller's dictionaries remain untouched. A supplied
    schema validates every run before any file is written and names the
    manifest's application. The manifest is written last; an explicit
    overwrite replaces each file atomically.
    """
    if isinstance(runs, (Mapping, str, bytes)):
        raise TypeError("runs must be an iterable of explicit run mappings")
    runs = list(runs)
    if not runs:
        raise ValueError("a job must contain at least one run")
    if baseseed is not None:
        if isinstance(baseseed, bool):
            raise TypeError("baseseed must be an integer")
        baseseed = operator.index(baseseed)
    prefix = Path(prefix)
    manifest = prefix.with_name(prefix.name + ".job.toml")
    documents = []
    entries = []
    resolved_runs = []
    for index, run in enumerate(runs):
        if not isinstance(run, Mapping):
            raise TypeError("each run must be a mapping of run sections")
        unknown = set(run) - set(_SECTIONS)
        if unknown:
            raise ValueError(f"unknown run sections: {sorted(unknown)}")
        sections = {name: _section(run.get(name)) for name in _SECTIONS}
        if baseseed is not None:
            sections["execution"].setdefault("seed", baseseed + index)
        document, resolved = _document(schema, sections)
        path = prefix.with_name(prefix.name + f".task{index + 1}.toml")
        documents.append((path, document))
        resolved_runs.append((path, resolved))
        entries.append("[[runs]]\nfile = " + _quoted(path.name) + "\n")
    header = "" if schema is None else "application = " + _quoted(resolved.application) + "\n\n"
    documents.append((manifest, header + "\n".join(entries)))
    _check_output_paths(resolved_runs, [path for path, _ in documents])
    _publish(documents, overwrite)
    return manifest


def read_job_manifest(filename):
    """Return the manifest's application (or None) and its absolute run files."""
    path = Path(filename)
    with open(path, "rb") as stream:
        document = tomllib.load(stream)
    unknown = set(document) - {"application", "runs"}
    if unknown:
        raise ValueError(f"unknown job-manifest keys: {sorted(unknown)}")
    application = document.get("application")
    if application is not None and (not isinstance(application, str) or not application):
        raise ValueError("job-manifest application must be a nonempty string")
    runs = document.get("runs")
    if not isinstance(runs, list) or not runs:
        raise ValueError("job manifest requires a nonempty [[runs]] array of run files")
    files = []
    for entry in runs:
        if not isinstance(entry, dict) or set(entry) != {"file"}:
            raise ValueError("each [[runs]] entry must contain only file")
        if not isinstance(entry["file"], str) or not entry["file"]:
            raise ValueError("[[runs]] file must be a nonempty filename")
        run = (path.parent / entry["file"]).resolve()
        if run in files:
            raise ValueError(f"duplicate TOML run file: {run}")
        files.append(run)
    return application, files


def _command(arguments, *, capture=False):
    from .tools import list2cmdline, log
    log(list2cmdline(arguments))
    return subprocess.run(arguments, check=True, text=True, capture_output=capture)


def execute(application, runs, *, mpi=None, mpirun="mpirun"):
    """Run TOML run files or job manifests with an application executable.

    Every run is validated by the application before the first one starts,
    and the runs then execute in order. Returns the absolute output.results
    path of each run. A failing process raises subprocess.CalledProcessError.
    """
    from .tools import check_existence
    application = os.fspath(application)
    runs = [Path(run).resolve() for run in ([runs] if isinstance(runs, (str, os.PathLike)) else runs)]
    if not runs:
        raise ValueError("at least one TOML run file is required")
    for path in runs:
        if path.suffix.lower() != ".toml":
            raise ValueError(f"expected a TOML run file or job manifest: {path}")
    check_existence(application)
    launcher = []
    if mpi is not None:
        if isinstance(mpi, bool) or operator.index(mpi) < 1:
            raise ValueError("mpi must be a positive process count")
        check_existence(mpirun)
        launcher = [os.fspath(mpirun), "-np", str(operator.index(mpi))]
    references, targets, seen = [], [], set()
    for path in runs:
        targets.append(path)
        expected, files = (read_job_manifest(path) if path.name.lower().endswith(".job.toml")
                           else (None, [path]))
        for file in files:
            if file.suffix.lower() != ".toml" or file.name.lower().endswith(".job.toml"):
                raise ValueError("job entries must name TOML run files")
            if file in seen:
                raise ValueError(f"duplicate TOML run file: {file}")
            seen.add(file)
            targets.append(file)
            schema = _command([application, "--schema", str(file)], capture=True).stdout
            run = run_config.load(file, schema)
            if expected is not None and expected != run.application:
                raise ValueError(f"job application {expected!r} does not match "
                                 f"executable application {run.application!r}")
            if "results" not in run.output:
                raise ValueError("the application schema must declare output.results")
            references.append((file, run))
    _check_output_paths(references, targets)
    # Scientific inputs are checked by the application before any run writes.
    for file, _ in references:
        _command([application, "--validate", str(file)])
    for file, _ in references:
        _command(launcher + [application, str(file)])
    return [str(Path(run.output["results"]).resolve()) for _, run in references]
