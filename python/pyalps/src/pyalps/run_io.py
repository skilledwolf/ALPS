# Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
"""Write schema-validated TOML runs and explicit multi-run job manifests.

The application supplies its schema separately. Arrays are parameter values;
each entry passed to ``write_run_files`` is one run, never an inferred sweep.
"""

from collections.abc import Mapping, Sequence
import json
import math
import operator
import os
from pathlib import Path
import tempfile

from . import run_config

_SECTIONS = ("parameters", "input", "output", "execution")


def _quoted(value):
    # TOML shares JSON's basic-string escapes but forbids surrogate escapes
    # and unescaped DEL. Keep complete Unicode characters in the UTF-8 file.
    return json.dumps(value, ensure_ascii=False).replace("\x7f", "\\u007f")


def _value(value):
    if isinstance(value, str):
        return _quoted(value)
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError("run values must be finite")
        return repr(value)
    if isinstance(value, complex):
        return "{ real = " + _value(value.real) + ", imag = " + _value(value.imag) + " }"
    # Native dictionary vectors return owning NumPy arrays or Python lists.
    if hasattr(value, "tolist"):
        value = value.tolist()
    if isinstance(value, Sequence):
        return "[" + ", ".join(_value(item) for item in value) + "]"
    raise TypeError(f"cannot serialize run value of type {type(value).__name__}")


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
    # Preserve relative paths in the document. The application's native loader
    # resolves them against this file, rather than the writer's working directory.
    resolved = run_config.resolve(schema, **supplied, base_directory="")
    lines = []
    for name in _SECTIONS:
        if supplied[name]:
            lines.append(f"[{name}]")
            values = getattr(resolved, name)
            for key in sorted(supplied[name]):
                lines.append(f"{_quoted(key)} = {_value(values[key])}")
            lines.append("")
    return "\n".join(lines), resolved


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
    """Publish complete files, refusing existing destinations by default."""
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
            descriptor, temporary = tempfile.mkstemp(dir=path.parent, prefix=".alps-run-")
            try:
                with os.fdopen(descriptor, "w", encoding="utf-8", newline="\n") as stream:
                    stream.write(document)
                if overwrite:
                    os.replace(temporary, path)
                else:
                    # Hard-link publication is atomic and cannot replace a file
                    # created by another writer since the initial existence check.
                    identity = Path(temporary).stat()
                    os.link(temporary, path)
                    created.append((path, identity))
            finally:
                Path(temporary).unlink(missing_ok=True)
    except BaseException:
        for path, identity in created:
            try:
                current = path.lstat()
                if (current.st_dev, current.st_ino) == (identity.st_dev, identity.st_ino):
                    path.unlink()
            except FileNotFoundError:
                pass
        raise


def write_run_file(filename, schema, *, parameters=None, input=None, output=None,
                   execution=None, overwrite=False):
    """Write one run and return its path, omitting schema-supplied defaults.

    ``schema`` is the application-owned TOML schema text. The same native
    validator used by C++ applications validates all four run sections.
    """
    path = Path(filename)
    document, resolved = _document(schema, dict(parameters=parameters, input=input,
                                                output=output, execution=execution))
    _check_output_paths([(path, resolved)], [path])
    _publish([(path, document)], overwrite)
    return path


def write_run_files(prefix, runs, schema, *, baseseed=None, overwrite=False):
    """Write explicit runs and return their ``<prefix>.job.toml`` manifest.

    Each run maps section names to mappings. If ``baseseed`` is supplied,
    missing ``execution.seed`` values receive ``baseseed + run_index``;
    explicit seeds and the caller's dictionaries remain untouched. The schema
    decides whether seeds are supported and which values are valid. All runs
    are validated before any file is written. Each file is published atomically;
    the manifest is published last. Explicit overwrite is atomic per file.
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
    documents.append((manifest, "application = " + _quoted(resolved.application) + "\n\n"
                      + "\n".join(entries)))
    _check_output_paths(resolved_runs, [path for path, _ in documents])
    _publish(documents, overwrite)
    return manifest
