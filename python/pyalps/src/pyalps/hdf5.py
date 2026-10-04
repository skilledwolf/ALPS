"""HDF5 IO through h5py, with explicit ownership transfer for native checkpoints.

Primitive values follow h5py/NumPy dtype and shape rules. Use NumPy arrays with
an explicit dtype for scientific data. Groups are ordinary h5py groups, without
list/dictionary guessing. Native callbacks receive NativeArchive instead; its
small primitive bridge accepts scalars and explicit NumPy arrays.
"""
# SPDX-License-Identifier: MIT
from contextlib import closing, contextmanager
import inspect
import os
import posixpath

import h5py

from .cxx.pyngshdf5_c import NativeArchive, register_archive_exception_type


class ArchiveError(Exception): pass
class ArchiveNotFound(ArchiveError, IOError): pass
class ArchiveClosed(ArchiveError, ValueError): pass
class InvalidPath(ArchiveError, SyntaxError): pass
class PathNotFound(ArchiveError, LookupError): pass
class WrongType(ArchiveError, TypeError): pass

for _index, _exception in enumerate((ArchiveError, ArchiveNotFound, ArchiveClosed,
                                   InvalidPath, PathNotFound, WrongType)):
    register_archive_exception_type(_index, _exception)
del register_archive_exception_type, _index, _exception


class archive:
    """Own one h5py file; native() transfers ownership for a complete operation.

    Modes are exactly r (read), a (create/update), and w (truncate). Native
    operations invalidate previously borrowed h5py groups/datasets, just as
    closing a h5py file does. Never retain the borrowed NativeArchive after the
    native() scope; copied callback views are closed when that scope exits.
    The two providers never share HDF5 identifiers.
    """
    def __init__(self, filename, mode="r"):
        if mode not in ("r", "a", "w"):
            raise ValueError("archive mode must be r, a, or w")
        self.filename = os.path.abspath(os.fspath(filename))
        self._mode = "r" if mode == "r" else "a"
        self._context = "/"
        self._native = None
        self._closed = False
        try:
            self._file = h5py.File(self.filename, mode)
        except FileNotFoundError as error:
            raise ArchiveNotFound(str(error)) from error

    @property
    def context(self):
        return self._context

    @property
    def is_open(self):
        return not self._closed

    @property
    def closed(self):
        return self._closed

    def __enter__(self):
        self._check_open()
        return self

    def __exit__(self, *exc):
        self.close()

    def __del__(self):
        try:
            self.close()
        except Exception:
            pass

    def _check_open(self):
        if self._closed:
            raise ArchiveClosed("I/O operation on closed file")

    def _owner(self):
        self._check_open()
        if self._native is not None:
            raise ArchiveError("native checkpoint IO currently owns this file")
        return self._file

    def close(self):
        if self._native is not None:
            raise ArchiveError("cannot close an archive inside native()")
        if not self._closed:
            self._closed = True
            self._file.close()

    def complete_path(self, path):
        path = os.fspath(path)
        return posixpath.normpath(path if path.startswith("/")
                                  else self.context + "/" + path)

    def set_context(self, path):
        self._owner()
        self._context = self.complete_path(path)

    def _location(self, path):
        path = self.complete_path(path)
        if "@" in path:
            obj, name = path.rsplit("@", 1)
            return obj.rstrip("/") or "/", name
        return path, None

    @contextmanager
    def native(self):
        """Borrow the native serializer's archive; h5py is closed until exit."""
        self._check_open()
        if self._native is not None:
            yield self._native
            return
        self._file.close()
        try:
            self._native = NativeArchive(self.filename, self._mode)
            self._native.set_context(self.context)
            yield self._native
        finally:
            try:
                if self._native is not None and self._native.is_open:
                    self._native.close()
            finally:
                self._native = None
                try:
                    self._file = h5py.File(self.filename, self._mode)
                except Exception:
                    self._closed = True
                    raise

    def __getitem__(self, path):
        owner = self._owner()
        obj, attr = self._location(path)
        try:
            node = owner[obj]
            if attr is not None:
                return node.attrs[attr]
            if isinstance(node, h5py.Group):
                return node
            if h5py.check_string_dtype(node.dtype) is not None:
                return node.asstr()[()]
            return node[()]
        except KeyError as error:
            raise PathNotFound(str(error)) from error

    def __setitem__(self, path, value):
        if getattr(value, "_alps_archive_savable", False) or inspect.ismethod(getattr(value, "save", None)):
            previous = self.context
            self.set_context(path)
            try:
                value.save(self)
            finally:
                self._context = previous
            return
        owner = self._owner()
        obj, attr = self._location(path)
        if hasattr(value, "save"):
            raise TypeError("object does not declare an archive-shaped save; call object.save(filename, path) instead")
        if attr is not None:
            owner[obj].attrs[attr] = value
            return
        # Create before replacing the old link: unsupported data cannot erase
        # the existing dataset or group. Dtypes are entirely h5py's concern.
        dataset = owner.create_dataset(None, data=value)
        parent, name = posixpath.split(obj)
        group = owner.require_group(parent or "/")
        if name in group:
            del group[name]
        group[name] = dataset

    def is_group(self, path):
        owner = self._owner()
        obj, attr = self._location(path)
        return attr is None and obj in owner and isinstance(owner[obj], h5py.Group)

    def is_data(self, path):
        owner = self._owner()
        obj, attr = self._location(path)
        return attr is None and obj in owner and isinstance(owner[obj], h5py.Dataset)

    def is_attribute(self, path):
        owner = self._owner()
        obj, attr = self._location(path)
        return attr is not None and obj in owner and attr in owner[obj].attrs

    def extent(self, path):
        owner = self._owner()
        obj, attr = self._location(path)
        if attr is None:
            shape = owner[obj].shape
        else:
            with closing(owner[obj].attrs.get_id(attr)) as identifier:
                shape = identifier.shape
        return None if shape is None else list(shape)

    def dimensions(self, path):
        shape = self.extent(path)
        return 0 if shape is None else len(shape)

    def is_scalar(self, path):
        return self.extent(path) == []

    def is_null(self, path):
        return self.extent(path) is None

    def is_complex(self, path):
        owner = self._owner()
        obj, attr = self._location(path)
        if attr is None:
            dtype = owner[obj].dtype
        else:
            with closing(owner[obj].attrs.get_id(attr)) as identifier:
                dtype = identifier.dtype
        return dtype.kind == "c"

    def list_children(self, path):
        return list(self._owner()[self.complete_path(path)])

    def list_attributes(self, path):
        obj, _ = self._location(path)
        return list(self._owner()[obj].attrs)

    def create_group(self, path):
        return self._owner().require_group(self.complete_path(path))

    def delete_data(self, path):
        del self._owner()[self.complete_path(path)]

    delete_group = delete_data

    def delete_attribute(self, path):
        obj, attr = self._location(path)
        if attr is None:
            raise InvalidPath("expected an attribute path")
        del self._owner()[obj].attrs[attr]
