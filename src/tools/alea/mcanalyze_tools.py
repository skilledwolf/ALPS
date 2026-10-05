# Copyright (C) 2011-2012 Lukas Gamper, Matthias Troyer, Maximilian Poprawe;
# 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
"""Shared command-line reader for native ALEA analysis results."""
from argparse import ArgumentParser
from pathlib import Path

import pyalps.alea as alea
import pyalps.hdf5 as h5
from pyalps import hdf5_name_decode, hdf5_name_encode


def impl_calculation(name, save_path, calculate):
    parser = ArgumentParser(description=f"Report native ALEA {name.lower()} estimates.")
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("-w", "--write", action="store_true", help="write estimates back to the archive")
    parser.add_argument("-n", "--name", action="append", dest="variables", metavar="VAR")
    parser.add_argument("-p", "--path", default="/simulation/results", metavar="HDF5-PATH")
    parser.add_argument("files", nargs="*")
    options = parser.parse_args()
    if not options.files:
        parser.print_help()
        return
    try:
        for filename in options.files:
            if not Path(filename).is_file():
                raise FileNotFoundError(filename)
            with h5.archive(filename, "a" if options.write else "r") as archive:
                variables = options.variables or [hdf5_name_decode(key) for key in archive.list_children(options.path)]
                estimates = []
                for variable in variables:
                    path = options.path.rstrip('/') + '/' + hdf5_name_encode(variable)
                    result = alea.read_result(archive, path)
                    estimates.append((variable, path, calculate(result)))
                # Validate all selected estimates before modifying this archive.
                for variable, path, estimate in estimates:
                    if options.verbose:
                        print(f"The {name} of variable {variable} in file {filename} is: {estimate}")
                    if options.write:
                        archive[path + '/' + save_path] = estimate
    except (ValueError, RuntimeError, LookupError, OSError, h5.ArchiveError) as error:
        parser.exit(1, f"{parser.prog}: {error}\n")
    print("Done")
