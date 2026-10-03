#!/usr/bin/env python3
"""Convert ALPS archive leaf encodings to ordinary HDF5 complex/Boolean types.

Requires h5py and NumPy, but no ALPS installation. Scientific group layouts are
preserved; this is not a solver checkpoint schema upgrade.
"""

# SPDX-License-Identifier: MIT

import argparse
from contextlib import closing
import itertools
import os
from pathlib import Path
import sys
import tempfile

import h5py
import numpy as np


COMPLEX = "__complex__"
TYPE = "__alps_type__"
BUFFER_BYTES = 8 * 1024 * 1024


def _attribute(obj, name):
    """Read bytes without h5py's automatic attribute string decoding."""
    with closing(obj.attrs.get_id(name)) as attr:
        if attr.shape is None:
            return h5py.Empty(attr.dtype)
        value = np.empty(attr.shape, attr.dtype)
        if value.dtype.hasobject:
            attr.read(value)
        else:
            with closing(attr.get_type()) as stored_type:
                attr.read(value, mtype=stored_type)
        return value


def _marker(name):
    return any(name == prefix or name.startswith(prefix + ":")
               for prefix in (COMPLEX, TYPE))


def _check_object(obj):
    """Reject structures whose meaning cannot survive an independent copy."""
    if isinstance(obj, h5py.Dataset):
        if obj.is_virtual or obj.external:
            raise ValueError(f"{obj.name}: virtual/external dataset storage is unsupported")
        with closing(obj.id.get_type()) as dtype:
            if dtype.detect_class(h5py.h5t.REFERENCE):
                raise ValueError(f"{obj.name}: object/region references are unsupported")
    for name in obj.attrs:
        with closing(obj.attrs.get_id(name)) as attr, closing(attr.get_type()) as dtype:
            if dtype.detect_class(h5py.h5t.REFERENCE):
                raise ValueError(f"{obj.name}@{name}: object/region references are unsupported")
        if _marker(name):
            if ":" in name:
                target = name.split(":", 1)[1]
                if target not in obj.attrs or _marker(target):
                    raise ValueError(f"{obj.name}@{name}: orphaned type marker")
            elif not isinstance(obj, h5py.Dataset):
                raise ValueError(f"{obj.name}@{name}: dataset marker on a group")


def _encoding(obj, attribute, boolean, report):
    suffix = "" if attribute is None else ":" + attribute
    where = obj.name if attribute is None else f"{obj.name}@{attribute}"
    if attribute is None:
        dtype, shape = obj.dtype, obj.shape
    else:
        with closing(obj.attrs.get_id(attribute)) as attr:
            dtype, shape = attr.dtype, attr.shape
    if shape is None:
        report.append(f"{where}: NULL dataspace preserved; original array rank is unavailable")
    kind = None
    if COMPLEX + suffix in obj.attrs:
        marker = _attribute(obj, COMPLEX + suffix)
        if (not isinstance(marker, np.ndarray) or marker.shape != ()
                or marker.dtype.kind not in "biu" or marker != 1):
            raise ValueError(f"{where}: complex marker must be scalar true")
        kind = "complex"
    if TYPE + suffix in obj.attrs:
        marker = _attribute(obj, TYPE + suffix)
        if not isinstance(marker, np.ndarray) or marker.shape != ():
            raise ValueError(f"{where}: type marker must be a scalar string")
        marker = marker.item()
        if isinstance(marker, bytes):
            marker = marker.decode("ascii")
        if marker not in ("bool", "int8") or kind is not None:
            raise ValueError(f"{where}: unknown or conflicting type markers")
        kind = marker
        if dtype.kind != "i" or dtype.itemsize != 1 or dtype.metadata:
            raise ValueError(f"{where}: {kind} marker requires plain signed-byte storage")
    if boolean:
        if kind not in (None, "bool"):
            raise ValueError(f"{where}: Boolean declaration conflicts with {kind} marker")
        if not (dtype.kind == "b" or (dtype.kind == "i" and dtype.itemsize == 1
                                      and not dtype.metadata)):
            raise ValueError(f"{where}: Boolean declaration requires signed-byte storage")
        kind = "bool"
    if kind == "complex":
        if dtype.kind not in "iuf" or dtype.metadata:
            raise ValueError(f"{where}: complex components must be plain real numbers")
        if shape is not None and (not shape or shape[-1] != 2):
            raise ValueError(f"{where}: complex payload must have a final dimension of two")
        dtype = np.dtype([("r", dtype), ("i", dtype)])
        shape = None if shape is None else shape[:-1]
    elif kind == "bool":
        dtype = np.dtype(bool)
    elif kind is None and dtype.kind == "i" and dtype.itemsize == 1 and not dtype.metadata:
        report.append(f"{where}: unmarked signed byte preserved as int8; "
                      "declare Boolean explicitly if appropriate")
    if kind is not None:
        report.append(f"{where}: {kind}")
    return kind, dtype, shape


def _values(data, kind, dtype, where):
    if isinstance(data, h5py.Empty):
        return h5py.Empty(dtype)
    data = np.asarray(data)
    if kind == "complex":
        value = np.empty(data.shape[:-1], dtype)
        # Arithmetic construction (real + 1j * imag) corrupts infinities and
        # signed zeros. Preserve each component, including integer precision.
        value["r"], value["i"] = data[..., 0], data[..., 1]
        return value
    if kind == "bool":
        if np.any((data != 0) & (data != 1)):
            raise ValueError(f"{where}: Boolean payload contains a value other than 0 or 1")
        return data.astype(bool)
    return data


def _blocks(shape, itemsize):
    """Tile every axis so even a very wide row does not require a huge buffer."""
    if not shape:
        yield ()
        return
    if 0 in shape:
        return
    budget = max(1, BUFFER_BYTES // itemsize)
    steps = []
    for extent in reversed(shape):
        step = min(extent, budget)
        steps.append(step)
        budget = max(1, budget // step)
    steps.reverse()
    for start in itertools.product(*(range(0, n, s) for n, s in zip(shape, steps))):
        yield tuple(slice(i, min(i + s, n)) for i, s, n in zip(start, steps, shape))


def _dataset(source, parent, name, kind, dtype, shape):
    if kind not in ("complex", "bool"):
        source.file.copy(source, parent, name=name, without_attrs=True)
        return parent[name]
    # Only lossless filters supported by h5py's public creation API are carried
    # across a datatype change. Never silently drop an unknown/lossy filter.
    with closing(source.id.get_create_plist()) as props:
        for index in range(props.get_nfilters()):
            if props.get_filter(index)[0] not in (1, 2, 3, 32000):
                raise ValueError(f"{source.name}: unsupported filter on converted dataset")
    options = {}
    if kind == "complex" and shape is not None and source.maxshape[-1] != 2:
        raise ValueError(f"{source.name}: extensible complex component axis is unsupported")
    if shape:
        if source.chunks is not None:
            options["chunks"] = source.chunks[:-1] if kind == "complex" else source.chunks
            options["maxshape"] = source.maxshape[:-1] if kind == "complex" else source.maxshape
            options.update(compression=source.compression, compression_opts=source.compression_opts,
                           shuffle=source.shuffle, fletcher32=source.fletcher32)
    if shape is not None:
        fill = source.fillvalue
        if kind == "complex":
            fill = np.full((2,), fill, dtype=source.dtype)
        options["fillvalue"] = _values(fill, kind, dtype, source.name)[()]
    target = parent.create_dataset(name, shape=shape, dtype=dtype, **options)
    if shape is not None:
        for selection in _blocks(shape, dtype.itemsize):
            source_selection = selection + (slice(None),) if kind == "complex" else selection
            target[selection] = _values(source[source_selection], kind, dtype, source.name)
    return target


def _address(obj):
    return h5py.h5o.get_info(obj.id).addr


def _copy(source, target, boolean_datasets, boolean_attributes):
    report = []
    declared = set()
    for path in boolean_datasets:
        obj = source[path]
        if not isinstance(obj, h5py.Dataset) or obj.file != source:
            raise ValueError(f"{path}: Boolean declaration must name a dataset in the input")
        declared.add((_address(obj), None))
    for path, name in boolean_attributes:
        obj = source[path]
        if name not in obj.attrs or _marker(name) or obj.file != source:
            raise ValueError(f"{path}@{name}: Boolean declaration must name an ordinary input attribute")
        declared.add((_address(obj), name))
    seen = {_address(source): target}

    def attributes(obj, out):
        for name in obj.attrs:
            if _marker(name):
                continue
            kind, dtype, shape = _encoding(obj, name, (_address(obj), name) in declared, report)
            value = _values(_attribute(obj, name), kind, dtype, f"{obj.name}@{name}")
            if kind in ("complex", "bool"):
                out.attrs.create(name, value, shape=shape, dtype=dtype)
            else:
                # NumPy dtype reconstruction loses HDF5 string padding/charset
                # and bitfield distinctions. Keep the original type and space.
                with closing(obj.attrs.get_id(name)) as attr, \
                        closing(attr.get_type()) as stored_type, \
                        closing(attr.get_space()) as space, \
                        closing(h5py.h5a.create(out.id, attr.get_name(), stored_type, space)) as copied:
                    if shape is not None:
                        copied.write(value, mtype=None if value.dtype.hasobject else stored_type)

    def group(obj, out):
        _check_object(obj)
        attributes(obj, out)
        for name in obj:
            link = obj.get(name, getlink=True)
            if isinstance(link, h5py.SoftLink):
                out[name] = link
                continue
            if not isinstance(link, h5py.HardLink):
                raise ValueError(f"{obj.name}/{name}: external links are unsupported")
            child = obj[name]
            address = _address(child)
            if address in seen:
                out[name] = seen[address]
                continue
            if isinstance(child, h5py.Group):
                result = out.create_group(name, track_order=True)
                seen[address] = result
                group(child, result)
            elif isinstance(child, h5py.Dataset):
                _check_object(child)
                kind, dtype, shape = _encoding(child, None, (address, None) in declared, report)
                result = _dataset(child, out, name, kind, dtype, shape)
                seen[address] = result
                attributes(child, result)
            else:
                raise ValueError(f"{child.name}: named datatypes are unsupported")

    group(source, target)
    return report


def convert(source, destination, *, boolean_datasets=(), boolean_attributes=()):
    """Write a new file; leave the source and any existing destination untouched."""
    source, destination = Path(source), Path(destination)
    if os.path.lexists(destination):
        raise FileExistsError(f"destination already exists: {destination}")
    temporary = None
    try:
        with h5py.File(source, "r") as src:
            if src.userblock_size:
                raise ValueError("files with a user block are unsupported")
            fd, name = tempfile.mkstemp(prefix=f".{destination.name}.", suffix=".tmp",
                                        dir=destination.absolute().parent)
            temporary = Path(name)
            os.close(fd)
            with h5py.File(temporary, "w", track_order=True) as dst:
                report = _copy(src, dst, boolean_datasets, boolean_attributes)
        # A sibling hard link publishes the complete file atomically and refuses
        # to overwrite a destination created by another process in the meantime.
        os.link(temporary, destination)
        return report
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    parser.add_argument("--boolean", action="append", default=[], metavar="DATASET",
                        help="interpret an unmarked signed-byte dataset as Boolean (repeatable)")
    parser.add_argument("--boolean-attribute", action="append", nargs=2, default=[],
                        metavar=("OBJECT", "ATTRIBUTE"),
                        help="interpret an unmarked signed-byte attribute as Boolean (repeatable)")
    args = parser.parse_args(argv)
    try:
        report = convert(args.source, args.destination, boolean_datasets=args.boolean,
                         boolean_attributes=args.boolean_attribute)
    except (OSError, ValueError, KeyError, TypeError) as error:
        parser.exit(1, f"{parser.prog}: {error}\n")
    for line in report:
        print(line)
    print(f"Converted {args.source} -> {args.destination}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
