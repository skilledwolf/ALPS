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
import threading
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
    protected = {path.resolve(): "a run or job file" for path in targets}
    paths = []
    for index, (path, resolved, schema) in enumerate(runs):
        # Without a schema, only the conventional result/checkpoint fields
        # are known paths. execute() obtains every run's native schema.
        rules = (tomllib.loads(schema) if schema is not None else
                 {"output": {key: {"type": "path"} for key in ("results", "checkpoint")}})
        for section in ("input", "output"):
            values = getattr(resolved, section)
            for key, rule in rules.get(section, {}).items():
                if rule.get("type") not in ("path", "path[]") or key not in values:
                    continue
                # Sidecar directories are used only when text output is enabled.
                if section == "output" and key == "text_directory" and (
                        "text" not in values or not values["text"]):
                    continue
                filenames = [values[key]] if rule["type"] == "path" else values[key]
                for filename in filenames:
                    destination = (path.parent / filename).resolve()
                    paths.append((destination, f"{section}.{key}", index))
                    if section == "input":
                        protected.setdefault(destination, f"input.{key} of {path.name}")
    owners = set()
    for destination, key, _ in paths:
        if not key.startswith("output."):
            continue
        if destination in protected:
            raise ValueError(f"{key} would overwrite {protected[destination]}: {destination}")
        if destination in owners:
            raise ValueError(f"multiple runs would write output file: {destination}")
        owners.add(destination)
    # Native validation checks each run's own sidecar names. Reserve active
    # text directories across runs without duplicating those application rules.
    for directory, key, owner in paths:
        if key != "output.text_directory":
            continue
        for destination, other_key, other_owner in paths:
            if other_owner != owner and destination.is_relative_to(directory):
                raise ValueError(f"{other_key} lies inside another run's "
                                 f"output.text_directory: {destination}")
    # Native snapshots own a clone/sweep filename namespace, even before the
    # first snapshot exists. Check it across separately executed job tasks.
    for prefix, key, _ in paths:
        if key != "output.snapshot_prefix":
            continue
        namespace = str(prefix) + ".clone"
        for destination in protected:
            if str(destination).startswith(namespace):
                raise ValueError(f"output.snapshot_prefix would overwrite {protected[destination]}")
        for destination, other_key, _ in paths:
            if str(destination).startswith(namespace):
                raise ValueError(f"output.snapshot_prefix overlaps {other_key}: {destination}")


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
    _check_output_paths([(path, resolved, schema)], [path])
    _publish([(path, document)], overwrite)
    return path


def write_run_files(prefix, runs, schema=None, *, baseseed=None, overwrite=False):
    """Write explicit runs and return their ``<prefix>.job.toml`` manifest.

    Each run maps section names to mappings. If ``baseseed`` is supplied,
    missing ``execution.seed`` values receive ``baseseed + run_index``;
    explicit seeds and the caller's dictionaries remain untouched. A supplied
    schema validates every run before any file is written and names the
    manifest's application. The manifest is written last; an explicit
    overwrite replaces each file atomically. Active text-output directories
    are reserved for their run; other runs' input and output paths must stay
    outside them.
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
        resolved_runs.append((path, resolved, schema))
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


def _run_all(commands, concurrency):
    """Run commands with at most ``concurrency`` at once, starting none after a failure."""
    if concurrency == 1:
        # An interrupt stops the active run, as subprocess.run does.
        for command in commands:
            _command(command)
        return
    pending, failures, lock = iter(commands), [], threading.Lock()

    def worker():
        while True:
            with lock:
                command = None if failures else next(pending, None)
            if command is None:
                return
            try:
                _command(command)
            except BaseException as error:
                with lock:
                    failures.append(error)

    workers = [threading.Thread(target=worker) for _ in range(min(concurrency, len(commands)))]
    for thread in workers:
        thread.start()
    try:
        for thread in workers:
            thread.join()
    except BaseException as error:
        # A terminal interrupt also reaches the runs in its process group.
        # Start no further run and wait for the active ones to stop.
        with lock:
            failures.insert(0, error)
        for thread in workers:
            thread.join()
    if failures:
        raise failures[0]


def execute(application, runs, *, mpi=None, mpirun="mpirun", concurrency=1):
    """Run TOML run files or job manifests with an application executable.

    Every run is validated by the application before the first one starts.
    Up to ``concurrency`` runs then execute at the same time, in order of
    submission, and their console output interleaves. Returns the absolute
    output.results path of each run in input order. A failing process raises
    subprocess.CalledProcessError after the other active runs finish; no
    further run starts. With ``concurrency > 1``, an interrupt also waits for
    the active runs. Active text-output directories cannot contain another
    run's input or output.

    With ``mpi``, each run starts its own ``mpirun -np mpi``. ``mpirun`` is an
    executable or an argument list. Concurrent MPI runs must not share cores:
    disable the launcher's default binding, for example with
    ``mpirun=["mpirun", "--bind-to", "none"]``.
    """
    from .tools import check_existence
    application = os.fspath(application)
    runs = [Path(run).resolve() for run in ([runs] if isinstance(runs, (str, os.PathLike)) else runs)]
    if not runs:
        raise ValueError("at least one TOML run file is required")
    for path in runs:
        if path.suffix.lower() != ".toml":
            raise ValueError(f"expected a TOML run file or job manifest: {path}")
    if isinstance(concurrency, bool) or operator.index(concurrency) < 1:
        raise ValueError("concurrency must be a positive number of runs")
    check_existence(application)
    launcher = []
    if mpi is not None:
        if isinstance(mpi, bool) or operator.index(mpi) < 1:
            raise ValueError("mpi must be a positive process count")
        launcher = ([os.fspath(mpirun)] if isinstance(mpirun, (str, os.PathLike))
                    else [os.fspath(argument) for argument in mpirun])
        if not launcher:
            raise ValueError("mpirun must name a launcher")
        check_existence(launcher[0])
        launcher += ["-np", str(operator.index(mpi))]
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
            references.append((file, run, schema))
    _check_output_paths(references, targets)
    # Scientific inputs are checked by the application before any run writes.
    for file, _, _ in references:
        _command([application, "--validate", str(file)])
    _run_all([launcher + [application, str(file)] for file, _, _ in references],
             operator.index(concurrency))
    return [str(Path(run.output["results"]).resolve()) for _, run, _ in references]
