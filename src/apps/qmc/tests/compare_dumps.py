#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compare optional QMC characterization dumps from two builds (requires h5py)."""
import argparse
from pathlib import Path

import h5py
import numpy as np


def datasets(archive):
    result = {}
    archive.visititems(
        lambda name, obj: result.update({name: obj})
        if isinstance(obj, h5py.Dataset) else None
    )
    return result


def compare(before, after):
    left = {p.name for p in before.glob("*.h5")}
    right = {p.name for p in after.glob("*.h5")}
    if not left or left != right:
        raise AssertionError(f"HDF5 file sets differ or are empty: {left ^ right}")
    count = 0
    for name in sorted(left):
        with h5py.File(before / name) as a, h5py.File(after / name) as b:
            ad, bd = datasets(a), datasets(b)
            if ad.keys() != bd.keys():
                raise AssertionError(f"{name}: dataset names differ")
            for path, av in ad.items():
                bv = bd[path]
                if av.shape != bv.shape or av.dtype != bv.dtype:
                    raise AssertionError(f"{name}:{path}: shape or dtype differs")
                aa, ba = np.asarray(av[()]), np.asarray(bv[()])
                # Object arrays hold variable-length HDF5 strings. Comparing
                # their pointer bytes would compare allocation addresses.
                equal = np.array_equal(aa, ba) if aa.dtype.hasobject else (
                    aa.tobytes() == ba.tobytes()
                )
                if not equal:
                    raise AssertionError(f"{name}:{path}: values differ")
                count += 1
    # Worm's preexisting uninitialized last_id_ also appears in kink IDs.
    # Its initialized solver state is compared in canonical HDF5 datasets.
    ax = {p.name for p in before.glob("*.xdr") if not p.name.startswith("worm_")}
    bx = {p.name for p in after.glob("*.xdr") if not p.name.startswith("worm_")}
    if ax != bx:
        raise AssertionError(f"XDR file sets differ: {ax ^ bx}")
    for name in sorted(ax):
        if (before / name).read_bytes() != (after / name).read_bytes():
            raise AssertionError(f"{name}: XDR bytes differ")
    print(f"Identical: {len(left)} HDF5 files, {count} datasets, {len(ax)} XDR payloads")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("before", type=Path)
    parser.add_argument("after", type=Path)
    args = parser.parse_args()
    compare(args.before, args.after)
