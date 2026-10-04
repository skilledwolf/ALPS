// Copyright (C) 2010-2012 Lukas Gamper; ALPS contributors.
// SPDX-License-Identifier: MIT
// Borrowed native checkpoint callbacks use this small scalar/NumPy bridge.
// Ordinary Python IO belongs to h5py; no Python container inference is done here.
#include <nanobind/nanobind.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/complex.hpp>
#include <alps/hdf5/pointer.hpp>
#include <alps/hdf5/vector.hpp>
#include "extract_from_pyobject.hpp"
#include "archive_savable.hpp"
#include "numpy_compat.hpp"
#include <algorithm>
#include <array>
#include <complex>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>
namespace nb = nanobind;
namespace alps::detail {
struct native_save_visitor {
    hdf5::archive & ar;
    std::string path;
    template<class T> void operator()(T const & value) const { ar[path] << value; }
    template<class T> void operator()(T const * values, std::vector<std::size_t> const & shape) const {
        if (shape.empty()) ar[path] << *values;
        else ar << make_pvp(path, values, shape);
    }
    void operator()(nb::list const &) const {
        throw nb::type_error("native callbacks require scalars or explicit NumPy arrays");
    }
    void operator()(nb::dict const &) const {
        throw nb::type_error("native callbacks require explicit scientific save/load methods");
    }
};
void native_save(hdf5::archive & ar, std::string const & relative, nb::handle data) {
    auto path=ar.complete_path(relative);
    if (nb::hasattr(data, "save") &&
        (nb::hasattr(data, pyalps::archive_savable_attr) ||
         alps::python::type_name(data.attr("save")) == "method")) {
        hdf5::detail::scoped_context context(ar, path);
        data.attr("save")(pyalps::owned_native_archive(ar));
    } else {
        if (nb::hasattr(data, "save"))
            throw nb::type_error("object does not declare an archive-shaped save; call object.save(filename, path) instead");
        // Unicode arrays have an explicit element family and rank; no
        // inference from Python lists or object arrays is needed.
        if (nb::isinstance(data, alps::python::numpy_module().attr("ndarray")) &&
            nb::cast<std::string>(data.attr("dtype").attr("kind")) == "U") {
            auto shape = nb::cast<std::vector<std::size_t>>(data.attr("shape"));
            auto values = nb::cast<std::vector<std::string>>(data.attr("ravel")().attr("tolist")());
            if (shape.empty()) ar.write(path, values.front());
            else ar.write(path, values.data(), shape);
        } else {
            native_save_visitor visitor{ar,path};
            extract_from_pyobject_py11(visitor,data);
        }
    }
}
    template <typename T>
    nb::object load_nd_array(alps::hdf5::archive & ar,
                             std::string const & path,
                             std::vector<std::size_t> const & shape) {
        std::size_t total = std::find(shape.begin(), shape.end(), 0) == shape.end() ? 1 : 0;
        if (total)
            for (auto size : shape) {
                if (total > std::numeric_limits<std::size_t>::max() / size)
                    throw hdf5::archive_error("HDF5 array extent overflows");
                total *= size;
            }
        nb::object array;
        if constexpr (std::is_same_v<T, bool>) {
            // vector<bool> has packed storage rather than a bool buffer.
            // Give NumPy ownership of the ordinary contiguous array read
            // by the archive, including a valid buffer for empty shapes.
            auto values = std::make_unique<bool[]>(total ? total : 1);
            if (shape.empty()) ar.read(path, values[0]);
            else if (total) ar.read(path, values.get(), shape);
            nb::capsule owner(values.get(), [](void * pointer) noexcept {
                delete[] static_cast<bool *>(pointer);
            });
            array = nb::cast(nb::ndarray<nb::numpy, bool>(
                values.release(), shape.size(), shape.data(), owner));
        } else {
            std::vector<T> flat(total);
            // Zero extents carry shape but have no payload to read.
            if (shape.empty()) ar.read(path, flat.front());
            else if (total) ar.read(path, flat.data(), shape);
            if constexpr (std::is_same_v<T, std::string>)
                array = alps::python::numpy_module().attr("array")(
                    nb::cast(flat), nb::arg("dtype") = "str").attr("reshape")(nb::cast(shape));
            else
                array = alps::python::make_numpy_array<T>(std::move(flat), shape);
        }
        return shape.empty() ? array.attr("__getitem__")(nb::tuple()) : array;
    }
    nb::object native_load(hdf5::archive & ar, std::string const & relative) {
        auto path=ar.complete_path(relative);
        if (ar.is_group(path))
            throw nb::type_error("read native checkpoint fields explicitly; groups have no Python container type");
        if (ar.is_null(path))
            throw hdf5::wrong_type(
                "NULL HDF5 value has no array shape; use the offline converter "
                "and an application schema to assign explicit empty dimensions");
        auto shape = ar.extent(path);
        if (shape.empty() && ar.is_datatype<std::string>(path)) {
            std::string value;
            ar[path] >> value;
            return nb::cast(value);
        }
        // One typed path preserves both NumPy scalar precision and array rank.
        #define TRY_TYPE(T) \
            if (ar.is_datatype<T>(path)) return load_nd_array<T>(ar, path, shape);
        TRY_TYPE(bool)
        TRY_TYPE(std::complex<float>)
        TRY_TYPE(std::complex<double>)
        TRY_TYPE(float)
        TRY_TYPE(double)
        TRY_TYPE(std::int8_t)
        TRY_TYPE(std::int16_t)
        TRY_TYPE(std::int32_t)
        TRY_TYPE(std::int64_t)
        TRY_TYPE(std::uint8_t)
        TRY_TYPE(std::uint16_t)
        TRY_TYPE(std::uint32_t)
        TRY_TYPE(std::uint64_t)
        TRY_TYPE(std::string)
        #undef TRY_TYPE
        throw hdf5::wrong_type("unsupported native callback datatype at " + path);
    }

    std::array<PyObject *, 6> exception_type = {};
    #define TRANSLATE_CPP_ERROR_TO_PYTHON(T, ID)                            \
        static void translate_ ## T (hdf5:: T const & e) {                  \
            std::string message =                                           \
                std::string(e.what()).substr(                               \
                    0, std::string(e.what()).find_first_of('\n'));          \
            PyErr_SetString(exception_type[ID] ? exception_type[ID]         \
                                               : PyExc_RuntimeError,       \
                            message.c_str());                               \
        }
    TRANSLATE_CPP_ERROR_TO_PYTHON(archive_error, 0)
    TRANSLATE_CPP_ERROR_TO_PYTHON(archive_not_found, 1)
    TRANSLATE_CPP_ERROR_TO_PYTHON(archive_closed, 2)
    TRANSLATE_CPP_ERROR_TO_PYTHON(invalid_path, 3)
    TRANSLATE_CPP_ERROR_TO_PYTHON(path_not_found, 4)
    TRANSLATE_CPP_ERROR_TO_PYTHON(wrong_type, 5)
    #undef TRANSLATE_CPP_ERROR_TO_PYTHON
    void register_exception_type(int id, nb::object type) {
        if (id < 0 || id >= static_cast<int>(exception_type.size()))
            throw std::out_of_range(
                "register_archive_exception_type: id out of range");
        // Keep a strong reference in the static table — the entry
        // is deliberately pinned until process exit because the
        // translators can fire at any time — but release any
        // previous entry so re-registration doesn't leak it.
        PyObject * previous = exception_type[id];
        Py_INCREF(type.ptr());
        exception_type[id] = type.ptr();
        Py_XDECREF(previous);
    }
}
NB_MODULE(pyngshdf5_c, m) {
    // Install the six C++→Python exception translators. Each calls the
    // matching translate_* above, which forwards to whichever Python
    // class was registered via register_archive_exception_type. If
    // pyalps/hdf5.py hasn't run yet, the translator falls back to
    // RuntimeError so the module is safely loadable on its own.
    nb::register_exception_translator(
        [](const std::exception_ptr &p, void * /*payload*/) {
            try { std::rethrow_exception(p); }
            catch (alps::hdf5::archive_not_found const & e) {
                alps::detail::translate_archive_not_found(e);
            } catch (alps::hdf5::archive_closed const & e) {
                alps::detail::translate_archive_closed(e);
            } catch (alps::hdf5::invalid_path const & e) {
                alps::detail::translate_invalid_path(e);
            } catch (alps::hdf5::path_not_found const & e) {
                alps::detail::translate_path_not_found(e);
            } catch (alps::hdf5::wrong_type const & e) {
                alps::detail::translate_wrong_type(e);
            } catch (alps::hdf5::archive_error const & e) {
                // Base class — must be caught LAST since the specialized
                // types above inherit from it.
                alps::detail::translate_archive_error(e);
            }
        });
    m.def("register_archive_exception_type",
          &alps::detail::register_exception_type);
    m.def("save_checkpoint", [](std::string const& filename, nb::object const& save) {
        alps::hdf5::save_checkpoint(filename, [&](auto& archive) {
            save(pyalps::owned_native_archive(archive));
        });
    }, nb::arg("filename"), nb::arg("save"));
    nb::class_<alps::hdf5::archive>(m, "NativeArchive")
        .def(nb::init<std::string, std::string>())
        .def_prop_ro("filename", &alps::hdf5::archive::get_filename)
        .def_prop_ro("context",  &alps::hdf5::archive::get_context)
        .def_prop_ro("is_open",  &alps::hdf5::archive::is_open)
        .def("set_context",     &alps::hdf5::archive::set_context)
        .def("is_group",        &alps::hdf5::archive::is_group)
        .def("is_data",         &alps::hdf5::archive::is_data)
        .def("is_attribute",    &alps::hdf5::archive::is_attribute)
        .def("close",           &alps::hdf5::archive::close)
        .def("extent",          &alps::hdf5::archive::extent)
        .def("dimensions",      &alps::hdf5::archive::dimensions)
        .def("is_scalar",       &alps::hdf5::archive::is_scalar)
        .def("is_complex",      &alps::hdf5::archive::is_complex)
        .def("is_null",         &alps::hdf5::archive::is_null)
        .def("list_children",   &alps::hdf5::archive::list_children)
        .def("list_attributes", &alps::hdf5::archive::list_attributes)
        .def("__setitem__",     &alps::detail::native_save)
        .def("__getitem__",     &alps::detail::native_load)
        .def("create_group",    &alps::hdf5::archive::create_group)
        .def("delete_data",     &alps::hdf5::archive::delete_data)
        .def("delete_group",    &alps::hdf5::archive::delete_group)
        .def("delete_attribute",&alps::hdf5::archive::delete_attribute);
}
