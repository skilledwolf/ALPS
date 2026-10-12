#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compare optional QMC characterization dumps from two builds (requires h5py)."""
import argparse
from pathlib import Path

import h5py
import numpy as np


def objects(archive):
    result = {"/": archive}
    archive.visititems(lambda name, obj: result.update({name: obj}))
    return result


def same_values(left, right):
    aa, ba = np.asarray(left), np.asarray(right)
    if aa.shape != ba.shape or aa.dtype != ba.dtype:
        return False
    # Object arrays hold variable-length HDF5 strings. Comparing their pointer
    # bytes would compare allocation addresses.
    return np.array_equal(aa, ba) if aa.dtype.hasobject else (
        aa.tobytes() == ba.tobytes()
    )


def compare(before, after):
    left = {p.name for p in before.glob("*.h5")}
    right = {p.name for p in after.glob("*.h5")}
    if not left or left != right:
        raise AssertionError(f"HDF5 file sets differ or are empty: {left ^ right}")
    count = 0
    attributes = 0
    for name in sorted(left):
        with h5py.File(before / name) as a, h5py.File(after / name) as b:
            ad, bd = objects(a), objects(b)
            if ad.keys() != bd.keys():
                raise AssertionError(f"{name}: group or dataset names differ")
            for path, av in ad.items():
                bv = bd[path]
                if isinstance(av, h5py.Dataset) != isinstance(bv, h5py.Dataset):
                    raise AssertionError(f"{name}:{path}: object types differ")
                if av.attrs.keys() != bv.attrs.keys():
                    raise AssertionError(f"{name}:{path}: attribute names differ")
                for attribute in av.attrs:
                    if not same_values(av.attrs[attribute], bv.attrs[attribute]):
                        raise AssertionError(f"{name}:{path}@{attribute}: values differ")
                    attributes += 1
                if not isinstance(av, h5py.Dataset):
                    continue
                if av.shape != bv.shape or av.dtype != bv.dtype:
                    raise AssertionError(f"{name}:{path}: shape or dtype differs")
                if not same_values(av[()], bv[()]):
                    raise AssertionError(f"{name}:{path}: values differ")
                count += 1
    # Worm's preexisting uninitialized last_id_ also appears in kink IDs.
    # Its initialized solver state is compared in canonical HDF5 datasets.
    ax = {p.name for p in before.glob("*.xdr")}
    bx = {p.name for p in after.glob("*.xdr")}
    if ax != bx:
        raise AssertionError(f"XDR file sets differ: {ax ^ bx}")
    deterministic = {name for name in ax if not name.startswith("worm_")}
    for name in sorted(deterministic):
        if (before / name).read_bytes() != (after / name).read_bytes():
            raise AssertionError(f"{name}: XDR bytes differ")
    print(f"Identical: {len(left)} HDF5 files, {count} datasets, "
          f"{attributes} attributes, {len(deterministic)} XDR payloads")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("before", type=Path)
    parser.add_argument("after", type=Path)
    args = parser.parse_args()
    compare(args.before, args.after)
