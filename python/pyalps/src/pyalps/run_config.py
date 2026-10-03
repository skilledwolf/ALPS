"""Application-owned TOML schemas with shared native validation.

``load`` resolves file paths relative to the run file. ``resolve`` keeps relative
paths unless ``base_directory`` is supplied. Arrays are parameter values; batches
of runs are managed separately by :mod:`pyalps.run_io`.
"""
# Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
import os
from .cxx.pyngsparams_c import (
    RunConfiguration,
    load_run_configuration,
    resolve_run_configuration,
)


def load(filename, schema):
    """Load a run using the schema supplied by its application."""
    return load_run_configuration(os.fspath(filename), schema)


def resolve(schema, *, parameters=None, input=None, output=None, execution=None,
            base_directory=None):
    """Validate owning snapshots of four section mappings and insert defaults."""
    return resolve_run_configuration(
        schema, dict(parameters or {}), dict(input or {}), dict(output or {}),
        dict(execution or {}), "" if base_directory is None else os.fspath(base_directory),
    )
