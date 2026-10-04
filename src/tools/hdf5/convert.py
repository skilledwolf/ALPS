#!/usr/bin/env python3
"""Convert ALPS archive leaf encodings to ordinary HDF5 complex/Boolean types.

Requires h5py and NumPy, but no ALPS installation. Scientific group layouts are
preserved. Typed alps.params.v1 dictionaries are upgraded to alps.params.v2;
other solver checkpoint schemas are not changed.
"""

# SPDX-License-Identifier: MIT

import argparse
from contextlib import closing
import itertools
import os
from pathlib import Path
import posixpath
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


def _copy_attribute(source, target, name):
    with closing(source.attrs.get_id(name)) as attr, closing(attr.get_type()) as datatype, \
            closing(attr.get_space()) as space, \
            closing(h5py.h5a.create(target.id, attr.get_name(), datatype, space)) as copied:
        if attr.shape is not None:
            value = _attribute(source, name)
            copied.write(value, mtype=None if value.dtype.hasobject else datatype)


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
        if data.dtype.kind == "b":
            data = data.view(np.uint8)
        if np.any((data != 0) & (data != 1)):
            raise ValueError(f"{where}: Boolean payload contains a value other than 0 or 1")
        return data.astype(bool)
    if kind == "params-text":
        return np.array([_text(value, where) for value in data.flat],
                        dtype=object).reshape(data.shape)
    return data


def _boolean_values(source, selection):
    if source.dtype.kind != "b":
        return source[selection]
    # h5py/NumPy's scalar Boolean conversion hides invalid raw enum codes.
    shape = tuple(part.stop - part.start for part in selection)
    values = np.empty(shape, dtype=np.uint8)
    with closing(source.id.get_type()) as datatype, closing(source.id.get_space()) as space:
        if shape:
            space.select_hyperslab(tuple(part.start for part in selection), shape)
        memory = h5py.h5s.create_simple(shape) if shape else h5py.h5s.create(h5py.h5s.SCALAR)
        with closing(memory):
            source.id.read(memory, space, values, mtype=datatype)
    return values


def _fill_value(source, kind, dtype):
    fill = source.fillvalue
    if kind == "bool" and source.dtype.kind == "b":
        fill = np.empty((), dtype=bool)
        with closing(source.id.get_create_plist()) as props:
            props.get_fill_value(fill)
    if kind == "complex":
        fill = np.full((2,), fill, dtype=source.dtype)
    return _values(fill, kind, dtype, source.name)[()]


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
    if kind == "params-text":
        dtype = h5py.string_dtype("utf-8")
    if source.shape is None and shape == (0,):
        # An explicit params vector type recovers the rank lost by NULL storage.
        return parent.create_dataset(name, shape=shape, dtype=dtype)
    if kind not in ("complex", "bool", "params-text"):
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
        options["fillvalue"] = _fill_value(source, kind, dtype)
    target = parent.create_dataset(name, shape=shape, dtype=dtype, **options)
    if shape is not None:
        for selection in _blocks(shape, BUFFER_BYTES if kind == "params-text" else dtype.itemsize):
            source_selection = selection + (slice(None),) if kind == "complex" else selection
            data = _boolean_values(source, source_selection) if kind == "bool" else source[source_selection]
            target[selection] = _values(data, kind, dtype, source.name)
    return target


def _address(obj):
    return h5py.h5o.get_info(obj.id).addr


def _schema_groups(source, pairs, matrices):
    groups = {}
    for kind, paths in (("pair", pairs), ("matrix", matrices)):
        for path in paths:
            path = posixpath.normpath("/" + str(path).lstrip("/"))
            obj, parents = source, set()
            for segment in path.strip("/").split("/") if path != "/" else ():
                if (not isinstance(obj, h5py.Group)
                        or not isinstance(obj.get(segment, getlink=True), h5py.HardLink)):
                    raise ValueError(f"{path}: schema selection must follow hard links")
                parents.add(_address(obj))
                obj = obj[segment]
            if not isinstance(obj, h5py.Group) or (kind == "matrix" and path == "/"):
                raise ValueError(f"{path}: {kind} selection must name a non-root group")
            address = _address(obj)
            if address in groups and groups[address][0] != kind:
                raise ValueError(f"{path}: conflicting schema selections")
            fields = {"first", "second"} if kind == "pair" else {
                "size1", "size2", "reserved_size1", "values"}
            if set(obj) != fields:
                raise ValueError(f"{path}: unexpected {kind} fields")
            if kind == "matrix":
                if h5py.h5o.get_info(obj.id).rc != 1:
                    raise ValueError(f"{path}: matrix group has hard-link aliases")
                for field in fields:
                    child = obj[field]
                    if (not isinstance(obj.get(field, getlink=True), h5py.HardLink)
                            or not isinstance(child, h5py.Dataset)
                            or h5py.h5o.get_info(child.id).rc != 1):
                        raise ValueError(f"{child.name}: matrix fields must be unaliased datasets")
            groups[address] = kind, path, parents
    selected = groups.keys()
    if any(parents & selected for _, _, parents in groups.values()):
        raise ValueError("overlapping schema selections are unsupported")
    return groups


def _matrix(group, parent):
    dimensions = []
    for field in ("size1", "size2", "reserved_size1"):
        dataset = group[field]
        with closing(dataset.id.get_type()) as datatype:
            integer = datatype.get_class() == h5py.h5t.INTEGER
        if (dataset.shape != () or dataset.dtype.kind not in "iu" or not integer
                or dataset.attrs):
            raise ValueError(f"{dataset.name}: matrix dimension must be an unannotated integer scalar")
        dimensions.append(int(dataset[()]))
    rows, columns, stride = dimensions
    values = group["values"]
    empty_null = values.shape is None and stride * columns == 0
    if (min(dimensions) < 0 or rows > stride
            or (values.shape != (stride * columns,) and not empty_null)):
        raise ValueError(f"{group.name}: invalid padded matrix dimensions")
    with closing(values.id.get_type()) as datatype:
        numeric_class = datatype.get_class() in (h5py.h5t.INTEGER, h5py.h5t.FLOAT,
                                                h5py.h5t.COMPOUND, h5py.h5t.ENUM)
    compound = (values.dtype.names == ("r", "i")
                and values.dtype.fields["r"][0] == values.dtype.fields["i"][0]
                and values.dtype.fields["r"][0].kind in "iuf")
    if not numeric_class or (values.dtype.kind not in "biufc" and not compound):
        raise ValueError(f"{values.name}: matrix values must be numeric")
    if values.maxshape != values.shape:
        raise ValueError(f"{values.name}: extensible matrix storage is unsupported")
    if set(group.attrs) & set(values.attrs):
        raise ValueError(f"{group.name}: matrix group/value attribute names conflict")
    with closing(values.id.get_create_plist()) as props:
        for index in range(props.get_nfilters()):
            if props.get_filter(index)[0] not in (1, 2, 3, 32000):
                raise ValueError(f"{values.name}: unsupported filter on converted matrix")
    boolean = values.dtype.kind == "b"
    options = {"fillvalue": _fill_value(values, "bool" if boolean else None, values.dtype)}
    if values.chunks is not None:
        options.update(chunks=True, compression=values.compression,
                       compression_opts=values.compression_opts, shuffle=values.shuffle,
                       fletcher32=values.fletcher32)
    target = parent.create_dataset(None, shape=(columns, rows), dtype=values.dtype, **options)
    for obj in (group, values):
        for name in obj.attrs:
            _copy_attribute(obj, target, name)
    for selection in _blocks(target.shape, values.dtype.itemsize):
        column_slice, row_slice = selection
        tile = np.empty((column_slice.stop - column_slice.start,
                         row_slice.stop - row_slice.start), dtype=values.dtype)
        for local, column in enumerate(range(column_slice.start, column_slice.stop)):
            flat_selection = (slice(column * stride + row_slice.start,
                                    column * stride + row_slice.stop),)
            data = (_boolean_values(values, flat_selection) if boolean else values[flat_selection])
            tile[local] = _values(data, "bool" if boolean else None, values.dtype, values.name)
        target[selection] = tile
    return target


def _text(value, where):
    value = value.decode("utf-8") if isinstance(value, bytes) else str(value)
    if "\0" in value:
        raise ValueError(f"{where}: string contains NUL")
    return value


def _scalar_text(dataset):
    if (not isinstance(dataset, h5py.Dataset) or dataset.shape != ()
            or h5py.check_string_dtype(dataset.dtype) is None):
        raise ValueError(f"{dataset.name}: expected a scalar string")
    return _text(dataset[()], dataset.name)


def _params_v1(source, declared, report):
    """Validate the explicit dictionary schema before publishing any output."""
    shapes, formats, payloads, texts = {}, set(), set(), set()

    def visit(_, obj):
        if not isinstance(obj, h5py.Group) or "format" not in obj:
            return
        format_dataset = obj["format"]
        if (not isinstance(format_dataset, h5py.Dataset) or format_dataset.shape != ()
                or h5py.check_string_dtype(format_dataset.dtype) is None):
            return
        if format_dataset[()] not in ("alps.params.v1", b"alps.params.v1"):
            return
        if "entries" not in obj or not isinstance(obj["entries"], h5py.Group):
            raise ValueError(f"{obj.name}: params checkpoint requires an entries group")
        entries = obj["entries"]
        if set(entries) != {str(i) for i in range(len(entries))}:
            raise ValueError(f"{entries.name}: params entry indices must be contiguous")
        names = set()
        for index in range(len(entries)):
            entry = entries[str(index)]
            if not isinstance(entry, h5py.Group) or set(entry) != {"name", "type", "value"}:
                raise ValueError(f"{entry.name}: invalid params entry")
            name, logical = _scalar_text(entry["name"]), _scalar_text(entry["type"])
            texts.update(_address(entry[field]) for field in ("name", "type"))
            payloads.update(_address(entry[field]) for field in ("name", "type", "value"))
            if name in names:
                raise ValueError(f"{entry.name}: duplicate parameter name")
            names.add(name)
            array = logical.endswith("[]")
            base = logical[:-2] if array else logical
            value = entry["value"]
            if not isinstance(value, h5py.Dataset):
                raise ValueError(f"{value.name}: params value must be a dataset")
            address = _address(value)
            boolean = base == "bool" or (address, None) in declared
            kind, dtype, shape = _encoding(value, None, boolean, [])
            with closing(value.id.get_type()) as datatype:
                physical_class = datatype.get_class()
            valid_type = {
                "bool": dtype.kind == "b",
                "int64": dtype.kind == "i" and dtype.itemsize == 8 and not dtype.metadata
                         and physical_class == h5py.h5t.INTEGER,
                "uint64": dtype.kind == "u" and dtype.itemsize == 8 and not dtype.metadata
                          and physical_class == h5py.h5t.INTEGER,
                "float64": dtype.kind == "f" and dtype.itemsize == 8 and not dtype.metadata
                           and physical_class == h5py.h5t.FLOAT,
                "complex128": (dtype.kind == "c" and dtype.itemsize == 16) or (
                    dtype.names == ("r", "i") and all(
                        dtype.fields[field][0].kind == "f" and dtype.fields[field][0].itemsize == 8
                        for field in ("r", "i"))),
                "string": h5py.check_string_dtype(dtype) is not None,
            }.get(base, False)
            if not valid_type:
                raise ValueError(f"{value.name}: params payload disagrees with declared {logical}")
            if array and shape is None:
                shapes[address] = (0,)
            elif shape is None or (len(shape) != 1 if array else shape != ()):
                raise ValueError(f"{value.name}: invalid params checkpoint rank")
            if base == "bool":
                declared.add((address, None))
            if base == "string":
                texts.add(address)
                if value.shape is not None:
                    for selection in _blocks(value.shape, BUFFER_BYTES):
                        for text in np.asarray(value[selection]).flat:
                            _text(text, value.name)
            # Validate Boolean values now as well as during streamed copying.
            if kind == "bool" and value.shape is not None:
                for selection in _blocks(value.shape, value.dtype.itemsize):
                    _values(_boolean_values(value, selection), kind, dtype, value.name)
        formats.add(_address(format_dataset))
        report.append(f"{obj.name}: alps.params.v1 -> alps.params.v2")

    visit("", source)
    source.visititems(visit)
    if formats & payloads:
        raise ValueError("params format dataset aliases an entry payload")
    return shapes, formats, texts


def _copy(source, target, boolean_datasets, boolean_attributes, pair_groups, matrix_groups):
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
    schemas = _schema_groups(source, pair_groups, matrix_groups)
    shapes, formats, texts = _params_v1(source, declared, report)
    seen = {_address(source): target}
    softlinks = []

    def attributes(obj, out):
        for name in obj.attrs:
            if _marker(name):
                continue
            kind, dtype, shape = _encoding(obj, name, (_address(obj), name) in declared, report)
            if kind in ("complex", "bool"):
                value = _values(_attribute(obj, name), kind, dtype, f"{obj.name}@{name}")
                out.attrs.create(name, value, shape=shape, dtype=dtype)
            else:
                # NumPy dtype reconstruction loses HDF5 string padding/charset
                # and bitfield distinctions. Keep the original type and space.
                _copy_attribute(obj, out, name)

    def group(obj, out):
        _check_object(obj)
        attributes(obj, out)
        for name in obj:
            link = obj.get(name, getlink=True)
            if isinstance(link, h5py.SoftLink):
                out[name] = link
                if schemas:
                    try:
                        original = obj[name]
                    except (KeyError, RuntimeError):
                        original = None
                    softlinks.append((_address(obj), name,
                                      _address(original) if original is not None else None))
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
                if address in shapes:
                    report.remove(f"{child.name}: NULL dataspace preserved; original array rank is unavailable")
                    report.append(f"{child.name}: declared empty vector shape (0,)")
                if address in formats:
                    result = out.create_dataset(name, data="alps.params.v2",
                                                dtype=h5py.string_dtype("utf-8"))
                else:
                    if address in texts:
                        kind = "params-text"
                    result = _dataset(child, out, name, kind, dtype, shapes.get(address, shape))
                seen[address] = result
                attributes(child, result)
            else:
                raise ValueError(f"{child.name}: named datatypes are unsupported")

    group(source, target)
    renamed = {}
    for address, (kind, path, _) in schemas.items():
        obj = seen[address]
        if kind == "pair":
            for old, new in (("first", "0"), ("second", "1")):
                obj.move(old, new)
                renamed[address, old] = new
        else:
            parent, name = target[posixpath.dirname(path)], posixpath.basename(path)
            null_storage = obj["values"].shape is None
            converted = _matrix(obj, parent)
            if null_storage:
                report.remove(f"{obj['values'].name}: NULL dataspace preserved; original array rank is unavailable")
                report.append(f"{path}: declared empty matrix shape {converted.shape}")
            del parent[name]
            parent[name] = converted
            seen[address] = converted
        report.append(f"{path}: legacy {kind} schema")
    for parent, name, original in softlinks:
        obj = seen[parent]
        name = renamed.get((parent, name), name)
        try:
            result = obj[name]
        except (KeyError, RuntimeError):
            result = None
        if ((original is None) != (result is None)
                or (original is not None and result.id != seen[original].id)):
            raise ValueError(f"{obj.name}/{name}: schema migration changes a soft-link target")
    return report


def convert(source, destination, *, boolean_datasets=(), boolean_attributes=(),
            pair_groups=(), matrix_groups=()):
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
                report = _copy(src, dst, boolean_datasets, boolean_attributes,
                               pair_groups, matrix_groups)
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
    parser.add_argument("--pair", action="append", default=[], metavar="GROUP",
                        help="migrate an explicitly selected first/second pair group (repeatable)")
    parser.add_argument("--matrix", action="append", default=[], metavar="GROUP",
                        help="migrate an explicitly selected padded numerical matrix group (repeatable)")
    args = parser.parse_args(argv)
    try:
        report = convert(args.source, args.destination, boolean_datasets=args.boolean,
                         boolean_attributes=args.boolean_attribute,
                         pair_groups=args.pair, matrix_groups=args.matrix)
    except (OSError, ValueError, KeyError, TypeError) as error:
        parser.exit(1, f"{parser.prog}: {error}\n")
    for line in report:
        print(line)
    print(f"Converted {args.source} -> {args.destination}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
