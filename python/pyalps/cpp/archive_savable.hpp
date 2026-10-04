// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
#ifndef PYALPS_ARCHIVE_SAVABLE_HPP
#define PYALPS_ARCHIVE_SAVABLE_HPP
#include <alps/hdf5/archive.hpp>
#include <nanobind/nanobind.h>
namespace pyalps {
inline constexpr char const * archive_savable_attr = "_alps_archive_savable";
inline void mark_archive_savable(nanobind::handle cls) {
    cls.attr(archive_savable_attr) = true;
}
// Native callbacks may retain their Python argument. Own a safe archive copy
// rather than exposing a reference to a serializer's stack object.
inline nanobind::object owned_native_archive(alps::hdf5::archive const & ar) {
    return nanobind::cast(new alps::hdf5::archive(ar), nanobind::rv_policy::take_ownership);
}
// Python h5py and native ALPS may use different HDF5 libraries. Transfer file
// ownership for the whole native operation, never a raw HDF5 object identifier.
// Virtual callbacks already borrowing a native archive stay with its owner.
template<class F> void with_native_archive(nanobind::handle object, F && operation) {
    namespace nb = nanobind;
    if (nb::isinstance<alps::hdf5::archive>(object)) {
        operation(nb::cast<alps::hdf5::archive &>(object));
        return;
    }
    if (!nb::hasattr(object, "native"))
        throw nb::type_error("expected an archive or NativeArchive");
    nb::object scope = object.attr("native")();
    nb::object native = scope.attr("__enter__")();
    try {
        operation(nb::cast<alps::hdf5::archive &>(native));
    } catch (...) {
        try { scope.attr("__exit__")(nb::none(), nb::none(), nb::none()); }
        catch (...) {} // preserve the scientific serializer's original error
        throw;
    }
    scope.attr("__exit__")(nb::none(), nb::none(), nb::none());
}
} // namespace pyalps
#endif
