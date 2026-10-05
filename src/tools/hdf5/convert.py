#!/usr/bin/env python3
"""Convert ALPS archive leaf encodings to ordinary HDF5 complex/Boolean types.

Requires h5py and NumPy, but no ALPS installation. Scientific group layouts are
preserved. Explicit profiles migrate released ALPS 3.0.0 parameter and ALEA
schemas; unspecified solver checkpoint schemas are not changed.
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
    if kind == "text":
        return np.array([_text(value, where) for value in data.flat],
                        dtype=object).reshape(data.shape)
    if kind == "numeric":
        return data.astype(dtype, copy=False)
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
    fill = _values(fill, kind, dtype, source.name)[()]
    if dtype.names == ("r", "i") and np.asarray(fill).dtype.names is None:
        # A NULL outer container's placeholder type may be integer. HDF5 does
        # not provide numeric-to-compound conversion for its default zero fill.
        if fill != 0:
            raise ValueError(f"{source.name}: NULL placeholder has a nonzero complex fill value")
        fill = np.zeros((), dtype=dtype)[()]
    return fill


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
    if kind == "text":
        dtype = h5py.string_dtype("utf-8")
    if source.shape is None and shape is not None:
        # A selected scientific schema recovers only known empty extents.
        return parent.create_dataset(name, shape=shape, dtype=dtype,
                                     fillvalue=_fill_value(source, kind, dtype))
    if kind not in ("complex", "bool", "text", "numeric"):
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
        for selection in _blocks(shape, BUFFER_BYTES if kind == "text" else dtype.itemsize):
            source_selection = selection + (slice(None),) if kind == "complex" else selection
            data = _boolean_values(source, source_selection) if kind == "bool" else source[source_selection]
            target[selection] = _values(data, kind, dtype, source.name)
    return target


def _address(obj):
    return h5py.h5o.get_info(obj.id).addr


def _schema_groups(source, pairs, matrices, parameters, alea, core_alea, alea_batches):
    groups = {}
    for kind, paths in (("pair", pairs), ("matrix", matrices),
                        ("parameters", parameters), ("alea", alea), ("alea-batches", alea_batches),
                        *(("core-alea:" + kind, [path]) for kind, path in core_alea)):
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
            if kind in ("pair", "matrix") and set(obj) != fields:
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
    value.encode("utf-8")
    if "\0" in value:
        raise ValueError(f"{where}: string contains NUL")
    return value


def _profiles(source, schemas, declared, report):
    """Normalize only types/extents established by a selected released schema."""
    conversions = {}

    def remember(dataset, kind, dtype, shape):
        address = _address(dataset)
        conversion = kind, np.dtype(dtype), shape
        if address in conversions and conversions[address] != conversion:
            raise ValueError(f"{dataset.name}: conflicting profile interpretations of a hard link")
        conversions[address] = conversion

    def encoding(dataset):
        _check_object(dataset)
        return _encoding(dataset, None, (_address(dataset), None) in declared, [])

    for kind, path, _ in schemas.values():
        group = source[path]
        if kind == "parameters":
            for name in group:
                dataset = group[name]
                if (not isinstance(group.get(name, getlink=True), h5py.HardLink)
                        or not isinstance(dataset, h5py.Dataset)):
                    raise ValueError(f"{group.name}/{name}: flat parameters require hard-linked datasets")
                _text(name, group.name)
                conversion, dtype, shape = encoding(dataset)
                with closing(dataset.id.get_type()) as datatype:
                    physical_class = datatype.get_class()
                if shape is not None and len(shape) > 1:
                    raise ValueError(f"{dataset.name}: parameters require scalars or rank-one vectors")
                if dtype.kind == "i" and dtype.itemsize == 1 and conversion not in ("bool", "int8"):
                    raise ValueError(f"{dataset.name}: ambiguous byte parameter; declare Boolean explicitly")
                if dtype.kind == "b":
                    conversion, dtype = "bool", np.dtype(bool)
                elif h5py.check_string_dtype(dtype) is not None:
                    conversion, dtype = "text", h5py.string_dtype("utf-8")
                elif dtype.kind in "iu" and physical_class == h5py.h5t.INTEGER and not dtype.metadata:
                    if dtype.itemsize > 8:
                        raise ValueError(f"{dataset.name}: unsupported parameter integer width")
                    conversion, dtype = "numeric", np.dtype("i8" if dtype.kind == "i" else "u8")
                elif dtype.kind == "f" and dtype.itemsize <= 8 and physical_class == h5py.h5t.FLOAT:
                    conversion, dtype = "numeric", np.dtype("f8")
                elif conversion == "complex" and dtype.fields["r"][0].kind == "f" and dtype.fields["r"][0].itemsize <= 8:
                    dtype = np.dtype([("r", "f8"), ("i", "f8")])
                elif dtype.kind == "c" and dtype.itemsize <= 16:
                    conversion, dtype = "numeric", np.dtype("c16")
                else:
                    raise ValueError(f"{dataset.name}: unsupported released parameter datatype")
                if shape is None:
                    shape = (0,)
                remember(dataset, conversion, dtype, shape)
            report.append(f"{path}: ALPS 3.0.0 flat parameters -> alps.params.v2")
        elif kind in ("alea", "alea-batches"):
            if "count" not in group or not isinstance(group["count"], h5py.Dataset):
                raise ValueError(f"{path}: ALEA profile requires a count dataset")
            count = group["count"]
            with closing(count.id.get_type()) as datatype:
                numeric = datatype.get_class() in (h5py.h5t.INTEGER, h5py.h5t.FLOAT)
            if (count.shape != () or not numeric or not np.isfinite(count[()])
                    or count[()] < 0 or count[()] != int(count[()])):
                raise ValueError(f"{count.name}: ALEA count must be a nonnegative integral scalar")
            for name in ("cannotrebin", "changed", "nonlinearoperations"):
                if name in group.attrs:
                    with closing(group.attrs.get_id(name)) as attribute:
                        if attribute.shape != ():
                            raise ValueError(f"{path}@{name}: ALEA flag must be scalar")
                    declared.add((_address(group), name))
            # Flagless mcdata and regular observable histories can look identical.
            # Do not guess the older sum-to-mean normalization convention.
            if ("mean/value" in group and "timeseries/data" in group
                    and not ("cannotrebin" in group.attrs or "changed" in group.attrs
                             or "timeseries/logbinning" in group)):
                raise ValueError(f"{path}: ambiguous pre-3.0 ALEA result flags/bin semantics")
            leaves = ("mean/value", "mean/error", "variance/value", "tau/value", "sum", "sum2",
                      "timeseries/partialbin", "timeseries/partialbin2")
            bins = ("timeseries/data", "timeseries/data2", "jacknife/data",
                    "timeseries/logbinning", "timeseries/logbinning2", "timeseries/logbinning_lastbin")
            element_shape, exemplars = None, {}
            for field in leaves + bins:
                if field not in group:
                    continue
                dataset = group[field]
                if not isinstance(dataset, h5py.Dataset):
                    raise ValueError(f"{dataset.name}: ALEA numeric field must be a dataset")
                conversion, dtype, shape = encoding(dataset)
                if shape is None:
                    continue
                candidate = shape[1:] if field in bins else shape
                if (len(candidate) > 1 or dtype.kind not in "iufc" and dtype.names != ("r", "i")
                        or dtype.kind in "iuf" and dtype.metadata):
                    raise ValueError(f"{dataset.name}: unsupported ALEA element shape/datatype")
                if element_shape is not None and element_shape != candidate:
                    raise ValueError(f"{dataset.name}: inconsistent ALEA element shapes")
                element_shape = candidate
                exemplars[field] = dtype
            for field in leaves + bins + ("mean/error_convergence", "timeseries/logbinning_counts", "labels"):
                if field not in group:
                    continue
                dataset = group[field]
                if not isinstance(dataset, h5py.Dataset):
                    raise ValueError(f"{dataset.name}: ALEA field must be a dataset")
                conversion, dtype, shape = encoding(dataset)
                if field == "labels":
                    if (h5py.check_string_dtype(dtype) is None
                            or shape is not None and (len(shape) > 1
                                or element_shape is not None and len(shape) != len(element_shape))):
                        raise ValueError(f"{dataset.name}: ALEA labels must match the scalar/vector element rank")
                    remember(dataset, "text", h5py.string_dtype("utf-8"), (0,) if shape is None else shape)
                    continue
                if shape is not None:
                    if field == "timeseries/logbinning_counts" and (len(shape) != 1 or dtype.kind not in "iu" or dtype.metadata):
                        raise ValueError(f"{dataset.name}: ALEA bin counts must be an integer vector")
                    if field == "mean/error_convergence" and (dtype.kind not in "iu" or dtype.metadata or shape != element_shape):
                        raise ValueError(f"{dataset.name}: invalid ALEA error convergence shape/type")
                    continue
                if field == "timeseries/logbinning_counts":
                    if dtype.kind not in "iu" or dtype.metadata:
                        raise ValueError(f"{dataset.name}: invalid ALEA bin counts datatype")
                    shape = (0,)
                elif element_shape is None:
                    raise ValueError(f"{dataset.name}: NULL ALEA storage lost its element shape/type")
                elif field in bins:
                    if element_shape:
                        # The outer empty vector stored INT NULL, losing T.
                        # mean is average_type<T>, which need not have T's dtype.
                        family = (("timeseries/data", "timeseries/data2", "timeseries/partialbin", "timeseries/partialbin2")
                                  if field in ("timeseries/data", "timeseries/data2") else
                                  ("mean/value", "sum", "timeseries/logbinning", "timeseries/logbinning2", "timeseries/logbinning_lastbin", "jacknife/data"))
                        dtype = next((exemplars[p] for p in family if p in exemplars), None)
                        if dtype is None:
                            raise ValueError(f"{dataset.name}: NULL ALEA bins lost their value datatype")
                    shape = (0,) + element_shape
                    conversion = None
                elif element_shape == (0,):
                    shape = (0,)
                    conversion = None
                else:
                    raise ValueError(f"{dataset.name}: NULL ALEA value disagrees with its scientific exemplar")
                remember(dataset, conversion, dtype, shape)
            report.append(f"{path}: ALPS 3.0.0 ALEA observable/result profile")
    return conversions


def _parameters(group, target):
    # Build anonymously so valid parameter names 'entries'/'format' never collide.
    entries = target.create_group(None, track_order=True)
    for name in list(group):
        entry = entries.create_group(str(len(entries)))
        entry.create_dataset("name", data=name, dtype=h5py.string_dtype("utf-8"))
        entry["value"] = group[name]
        del group[name]
    group["entries"] = entries
    group.create_dataset("format", data="alps.params.v2", dtype=h5py.string_dtype("utf-8"))


CORE_ALEA_KINDS = {"mean": 1, "variance": 2, "covariance": 3, "autocorr": 4, "batch": 5}


def _core_alea(group, family):
    """Tag explicitly selected released result layouts; never invent missing state."""
    if family not in CORE_ALEA_KINDS:
        raise ValueError(f"{group.name}: unknown Core ALEA result kind {family}")
    if "version" in group.attrs or "kind" in group.attrs:
        raise ValueError(f"{group.name}: expected an unversioned released Core ALEA result")

    def uint64(identifier):
        with closing(identifier.get_type()) as datatype:
            return (datatype.get_class() == h5py.h5t.INTEGER and datatype.get_size() == 8
                    and datatype.get_sign() == h5py.h5t.SGN_NONE)

    def integer(name, attribute=False):
        value = _attribute(group, name) if attribute else np.asarray(group[name][()])
        if attribute:
            with closing(group.attrs.get_id(name)) as identifier:
                canonical = uint64(identifier)
        else: canonical = uint64(group[name].id)
        if value.shape != () or not canonical:
            raise ValueError(f"{group.name}/{name}: expected a scalar uint64")
        return int(value)

    size = integer("size", True)
    if not size:
        raise ValueError(f"{group.name}: zero-component ALEA results are unsupported")

    def shaped(name, shapes, complex_allowed=True):
        value = group[name]
        if (not isinstance(value, h5py.Dataset) or value.shape not in shapes
                or not (value.dtype.kind == "f" and value.dtype.itemsize == 8
                        or complex_allowed and value.dtype.kind == "c" and value.dtype.itemsize == 16)):
            raise ValueError(f"{group.name}/{name}: invalid {family} result shape/datatype")

    shaped("mean/value", [(size,)])
    if family != "mean": shaped("mean/error", [(size,), (size, 2, 2)], False)
    if family == "autocorr":
        levels = integer("nlevel", True)
        if not levels or set(group["level"]) != {str(i) for i in range(levels)}:
            raise ValueError(f"{group.name}: invalid autocorrelation levels")
        for i in range(levels):
            level = group[f"level/{i}"]
            if _attribute(level, "size") != size:
                raise ValueError(f"{level.name}: inconsistent autocorrelation component count")
            _core_alea(level, "variance")
    elif family == "batch":
        batches = integer("num_batches", True)
        shaped("batch/sum", [(batches, size)])
        counts = group["batch/count"]
        if not batches or counts.shape != (batches,) or not uint64(counts.id):
            raise ValueError(f"{group.name}: invalid per-batch counts")
    else:
        count = integer("count")
        if family in ("variance", "covariance"):
            count2 = group["count2"]
            if count2.shape != () or count2.dtype.kind != "f" or count2.dtype.itemsize != 8:
                raise ValueError(f"{group.name}: invalid squared-count statistic")
            value = count2[()]
            if not np.isfinite(value) or value < 0 or (count == 0) != (value == 0):
                raise ValueError(f"{group.name}: inconsistent squared-count statistic")
            if family == "variance": shaped("var", [(size,), (size, 2, 2)], False)
            else: shaped("cov", [(size, size), (size, size, 2, 2)])
    group.attrs["version"] = np.uint64(1)
    group.attrs["kind"] = np.uint32(CORE_ALEA_KINDS[family])


def _alea_batches(group):
    """Recover released linear bin sums/counts, never an accumulator cursor."""
    if "version" in group.attrs or "kind" in group.attrs or "batch" in group:
        raise ValueError(f"{group.name}: expected a released ALPS 3.0.0 ALEA group")

    def integer(value, where):
        value = np.asarray(value)
        if (value.shape != () or value.dtype.kind not in "iuf" or not np.isfinite(value)
                or value < 0 or value != int(value) or int(value) > np.iinfo("u8").max):
            raise ValueError(f"{where}: expected a nonnegative uint64 count")
        return int(value)

    count = integer(group["count"][()], group.name + "/count")
    # DetailedBinning and SimpleObservableData store sums; MCData stores means.
    raw = "timeseries/logbinning" in group
    mcdata = "cannotrebin" in group.attrs
    evaluator = "changed" in group.attrs and "nonlinearoperations" in group.attrs
    if sum((raw, mcdata, evaluator)) != 1:
        raise ValueError(f"{group.name}: ambiguous or summary-only ALEA bin semantics")
    if any(bool(group.attrs.get(flag, False))
           for flag in ("cannotrebin", "changed", "nonlinearoperations")):
        raise ValueError(f"{group.name}: nonlinear/transformed ALEA bins cannot be recovered")
    if "timeseries/data" not in group:
        raise ValueError(f"{group.name}: missing linear bin history")
    data = group["timeseries/data"]
    if (data.shape is None or len(data.shape) not in (1, 2)
            or data.dtype.kind not in "fc" or data.dtype.itemsize > 16
            or _text(data.attrs.get("binningtype", ""), data.name) != "linear"):
        raise ValueError(f"{data.name}: expected scalar/vector floating-point linear bins")
    size = data.shape[1] if len(data.shape) == 2 else 1
    if not size:
        raise ValueError(f"{group.name}: zero-component bins are unsupported")
    binsize = integer(data.attrs.get("binsize", -1), data.name + "@binsize")
    discarded = integer(data.attrs.get("discard", 0), data.name + "@discard")
    if discarded:
        raise ValueError(f"{data.name}: discarded bins cannot reconstruct the full sample count")
    full = data.shape[0]
    partial, remainder = None, 0
    if raw and "timeseries/partialbin" in group:
        partial = group["timeseries/partialbin"]
        if partial.shape != data.shape[1:] or partial.dtype != data.dtype:
            raise ValueError(f"{partial.name}: inconsistent partial-bin shape/datatype")
        remainder = integer(partial.attrs.get("count", -1), partial.name + "@count")
    # Evaluator/MCData writers include the final unfinished row in data itself;
    # its weight is fixed by the recorded total count and preceding full bins.
    last = count - (full - 1) * binsize if full and not raw else binsize
    complete = count == full * binsize + remainder if raw else (
        0 < last <= binsize if full else count == 0)
    if not complete or (count or full) and not binsize or remainder > binsize:
        raise ValueError(f"{group.name}: incomplete/inconsistent bin counts; missing samples cannot be invented")

    # Owned legacy leaves change shape or disappear. Reject aliases to those
    # objects so conversion cannot silently break a hard-link graph.
    leaves = ("count", "sum", "sum2", "mean/value", "mean/error", "mean/error_convergence",
              "variance/value", "tau/value", "jacknife/data", "timeseries/data",
              "timeseries/data2", "timeseries/partialbin", "timeseries/partialbin2",
              "timeseries/logbinning", "timeseries/logbinning2",
              "timeseries/logbinning_lastbin", "timeseries/logbinning_counts")
    for path in leaves:
        if path not in group:
            continue
        current = group
        for segment in path.split("/"):
            if not isinstance(current.get(segment, getlink=True), h5py.HardLink):
                raise ValueError(f"{group.name}/{path}: ALEA state must follow hard links")
            current = current[segment]
            if h5py.h5o.get_info(current.id).rc != 1:
                raise ValueError(f"{current.name}: ALEA state has hard-link aliases")

    batches = max(2, full + bool(remainder))
    dtype = np.dtype("c16" if data.dtype.kind == "c" else "f8")
    counts = group.create_dataset("batch/count", shape=(batches,), dtype="u8")
    sums = group.create_dataset("batch/sum", shape=(batches, size), dtype=dtype)
    counts[:full] = binsize
    if full and not raw:
        counts[full - 1] = last
    if remainder:
        counts[full] = remainder
    # Tile both dimensions; this also bounds memory for unusually wide vectors.
    for selection in _blocks(data.shape, dtype.itemsize):
        values = data[selection].astype(dtype)
        if mcdata:
            values *= binsize
        if not np.isfinite(values).all():
            raise ValueError(f"{data.name}: bin sums must be finite")
        target = selection if len(selection) == 2 else selection + (slice(0, 1),)
        sums[target] = values if len(selection) == 2 else values[:, None]
    if remainder:
        for selection in _blocks(partial.shape, dtype.itemsize):
            values = np.asarray(partial[selection], dtype=dtype)
            if not np.isfinite(values).all():
                raise ValueError(f"{partial.name}: partial-bin sums must be finite")
            sums[(full,) + selection if selection else (full, 0)] = values

    # Use the native weighted-bin estimator. Published legacy error estimates
    # may differ; raw bins preserve the information for native joint analysis.
    mean = group.create_dataset("batch/mean", shape=(size,), dtype=dtype)
    error = group.create_dataset("batch/error", shape=(size,), dtype="f8")
    count2 = full * binsize**2 + remainder**2 if raw else (
        (full - 1) * binsize**2 + last**2 if full else 0)
    for (components,) in _blocks((size,), dtype.itemsize * 3):
        total = np.zeros(components.stop - components.start, dtype=dtype)
        for rows, _ in _blocks((batches, len(total)), dtype.itemsize):
            total += sums[rows, components].sum(axis=0)
        average = total / count if count else np.full(total.shape, np.nan, dtype=dtype)
        if raw and "sum" in group:
            recorded = group["sum"]
            if recorded.shape != data.shape[1:]:
                raise ValueError(f"{recorded.name}: inconsistent total-sum shape")
            expected = np.asarray(recorded[components] if recorded.shape else recorded[()]).reshape(-1)
            precision = max(1e-12, 32 * np.finfo(data.dtype).eps)
            if not np.isfinite(expected).all() or not np.allclose(total, expected, rtol=precision, atol=precision):
                raise ValueError(f"{group.name}: bin sums disagree with the recorded total sum")
        variance = np.zeros(total.shape, dtype="f8")
        for rows, _ in _blocks((batches, len(total)), dtype.itemsize * 3):
            weight = counts[rows]
            occupied = weight != 0
            values = sums[rows, components][occupied] / weight[occupied, None]
            variance += (weight[occupied, None] * abs(values - average)**2).sum(axis=0)
        if not count:
            uncertainty = np.full(total.shape, np.nan)
        elif count**2 == count2:
            uncertainty = np.full(total.shape, np.inf)
        else:
            uncertainty = np.sqrt(variance * count2 / (count * (count**2 - count2)))
        mean[components], error[components] = average, uncertainty
    for path in leaves:
        if path in group:
            del group[path]
    for name in ("variance", "tau", "jacknife", "timeseries", "mean"):
        if name in group and not len(group[name]) and not len(group[name].attrs):
            del group[name]
    group.require_group("mean")
    group.move("batch/mean", "mean/value")
    group.move("batch/error", "mean/error")
    for flag in ("cannotrebin", "changed", "nonlinearoperations"):
        if flag in group.attrs:
            del group.attrs[flag]
    group.attrs["size"], group.attrs["num_batches"] = np.uint64(size), np.uint64(batches)
    _core_alea(group, "batch")


def _copy(source, target, boolean_datasets, boolean_attributes, pair_groups, matrix_groups,
          parameter_groups, alea_groups, core_alea_groups, alea_batch_groups):
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
    schemas = _schema_groups(source, pair_groups, matrix_groups, parameter_groups, alea_groups,
                             core_alea_groups, alea_batch_groups)
    conversions = _profiles(source, schemas, declared, report)
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
                if address in conversions and shape is None:
                    report.remove(f"{child.name}: NULL dataspace preserved; original array rank is unavailable")
                    report.append(f"{child.name}: schema establishes empty shape {conversions[address][2]}")
                kind, dtype, shape = conversions.get(address, (kind, dtype, shape))
                result = _dataset(child, out, name, kind, dtype, shape)
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
        elif kind == "matrix":
            parent, name = target[posixpath.dirname(path)], posixpath.basename(path)
            null_storage = obj["values"].shape is None
            converted = _matrix(obj, parent)
            if null_storage:
                report.remove(f"{obj['values'].name}: NULL dataspace preserved; original array rank is unavailable")
                report.append(f"{path}: declared empty matrix shape {converted.shape}")
            del parent[name]
            parent[name] = converted
            seen[address] = converted
        elif kind == "parameters":
            _parameters(obj, target)
        elif kind == "alea-batches":
            _alea_batches(obj)
            report.append(f"{path}: ALPS 3.0.0 linear bins -> native ALEA batch result; error recomputed")
        elif kind.startswith("core-alea:"):
            _core_alea(obj, kind.split(":", 1)[1])
            report.append(f"{path}: ALPSCore 2.3.3 ALEA result -> versioned native result")
        if kind in ("pair", "matrix"):
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
            pair_groups=(), matrix_groups=(), parameter_groups=(), alea_groups=(), core_alea_groups=(),
            alea_batch_groups=()):
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
                               pair_groups, matrix_groups, parameter_groups, alea_groups, core_alea_groups,
                               alea_batch_groups)
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
    parser.add_argument("--parameters", action="append", default=[], metavar="GROUP",
                        help="migrate ALPS 3.0.0 flat parameters to native typed checkpoints (repeatable)")
    parser.add_argument("--alea", action="append", default=[], metavar="GROUP",
                        help="normalize one ALPS 3.0.0 ALEA observable/result (repeatable)")
    parser.add_argument("--alea-batches", action="append", default=[], metavar="GROUP",
                        help="recover complete ALPS 3.0.0 linear bins as a native ALEA analysis result; "
                             "recompute uncertainty, never invent restart state (repeatable)")
    parser.add_argument("--core-alea", action="append", nargs=2, default=[],
                        metavar=("KIND", "GROUP"),
                        help="migrate a released ALPSCore 2.3.3 ALEA result; KIND is "
                             + ", ".join(CORE_ALEA_KINDS) + " (repeatable)")
    args = parser.parse_args(argv)
    try:
        report = convert(args.source, args.destination, boolean_datasets=args.boolean,
                         boolean_attributes=args.boolean_attribute,
                         pair_groups=args.pair, matrix_groups=args.matrix,
                         parameter_groups=args.parameters, alea_groups=args.alea,
                         core_alea_groups=args.core_alea, alea_batch_groups=args.alea_batches)
    except (OSError, ValueError, KeyError, TypeError) as error:
        parser.exit(1, f"{parser.prog}: {error}\n")
    for line in report:
        print(line)
    print(f"Converted {args.source} -> {args.destination}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
