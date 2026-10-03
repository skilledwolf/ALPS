// Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>
//               Matthias Troyer <troyer@comp-phys.org>
//               2026       by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
#include <nanobind/nanobind.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/vector.h>
#include <alps/hdf5/archive.hpp>
#include "archive_savable.hpp"
#include <alps/params.hpp>
#include <complex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <boost/variant/apply_visitor.hpp>
#include <boost/variant/static_visitor.hpp>
#include <vector>
#include "dict_to_params.hpp"
#include "run_config.hpp"
namespace nb = nanobind;
namespace {
// Walk the paramvalue variant and wrap each native alternative as a
// nb::object. Called from __getitem__.
struct paramvalue_to_py_visitor : boost::static_visitor<nb::object> {
    nb::object operator()(alps::params_ns::detail::None) const { return nb::none(); }
    template <typename T>
    nb::object operator()(T const & value) const {
        return nb::cast(value);
    }
    template <typename T>
    nb::object operator()(std::vector<T> const & value) const {
        // String parameters need variable-length elements: NumPy's inferred
        // fixed-width Unicode dtype silently truncates longer replacements.
        if constexpr (std::is_same<T, std::string>::value)
            return nb::cast(value);
        // An empty numeric sequence has no elements from which NumPy can
        // infer its type. Preserve the native family, including Boolean masks.
        else
            return alps::python::numpy_module().attr("array")(
                nb::cast(value), nb::arg("dtype") = alps::python::numpy_dtype<T>::name);
    }
};
nb::object paramvalue_to_py(const alps::params::value_type& value) {
    return value.apply_visitor(paramvalue_to_py_visitor());
}
// Assigning copies the supported value; callers reassign after editing arrays.
void params_setitem(alps::params & self, nb::object const & key_obj, nb::handle value) {
    pyalps::set_param_value(self, nb::cast<std::string>(nb::str(key_obj)), value);
}
nb::object params_getitem(alps::params & self, nb::object const & key_obj) {
    std::string key = nb::cast<std::string>(nb::str(key_obj));
    const auto it=self.find(key);
    if(it==self.end() || it->second.empty()) throw nb::key_error(key.c_str());
    return paramvalue_to_py(it->second);
}
void params_delitem(alps::params & self, nb::object const & key_obj) {
    const auto key=nb::cast<std::string>(nb::str(key_obj));
    if(!self.exists(key)) throw nb::key_error(key.c_str());
    self.erase(key);
}
bool params_contains(alps::params & self, nb::object const & key_obj) {
    return self.exists(nb::cast<std::string>(nb::str(key_obj)));
}
nb::object value_or_default(alps::params & self, nb::object const & key, nb::object const & dflt) {
    return params_contains(self, key) ? params_getitem(self, key) : dflt;
}
void params_load(alps::params & self, alps::hdf5::archive & ar, std::string const & path) {
    alps::hdf5::archive reader(ar);
    reader.set_context(ar.complete_path(path));
    self.load(reader);
}
std::string params_print(alps::params & self) {
    std::stringstream ss;
    ss << self;
    return ss.str();
}
// All values own native storage; a dictionary copy is already deep.
alps::params params_deepcopy(const alps::params& self, nb::handle) { return self; }
}  // namespace
NB_MODULE(pyngsparams_c, m) {
    nb::class_<alps::params>(m, "params")
        .def("__init__", [](alps::params * self) {
            new (self) alps::params(pyalps::params_from_dict(nb::dict()));
        })
        .def("__init__",
             [](alps::params * self, nb::dict const & d) {
                 new (self) alps::params(pyalps::params_from_dict(d));
             },
             nb::arg("dict"))
        .def("__init__", [](alps::params * self, alps::hdf5::archive & ar, std::string const & path) {
                 alps::params loaded;
                 params_load(loaded, ar, path);
                 new (self) alps::params(loaded);
             },
             nb::arg("archive"),
             nb::arg("path") = std::string("/parameters"))
        .def("__eq__", [](const alps::params& a, const alps::params& b) { return a == b; }, nb::is_operator())
        .def("__ne__", [](const alps::params& a, const alps::params& b) { return a != b; }, nb::is_operator())
        .def("__eq__", [](const alps::params& a, const nb::dict& b) { return a == pyalps::params_from_dict(b); }, nb::is_operator())
        .def("__ne__", [](const alps::params& a, const nb::dict& b) { return a != pyalps::params_from_dict(b); }, nb::is_operator())
        .def("__len__",      [](alps::params const & self) { return self.size(); })
        .def("__deepcopy__", &params_deepcopy)
        .def("__getitem__",  &params_getitem)
        .def("__setitem__",  &params_setitem, nb::arg("key"), nb::arg("value").none())
        .def("__delitem__",  &params_delitem)
        .def("__contains__", &params_contains)
        .def("__iter__",     [](alps::params & self) {
                                 // Snapshot keys: params' native iterator
                                 // can be invalidated by deletion even while the
                                 // params object itself remains alive.
                                 nb::list keys;
                                 for (auto const & entry : self)
                                     keys.append(nb::cast(entry.first));
                                 return keys.attr("__iter__")();
                             })
        .def("__str__",      &params_print)
        .def("valueOrDefault", &value_or_default)
        .def("save",         &alps::params::save)
        .def("load",         &params_load,
             nb::arg("archive"),
             nb::arg("path") = std::string("/parameters"));
    pyalps::mark_archive_savable(m.attr("params"));
    nb::class_<alps::run_configuration>(m, "RunConfiguration")
        .def_ro("application", &alps::run_configuration::application)
        .def_ro("schema_version", &alps::run_configuration::schema_version)
        .def_ro("parameters", &alps::run_configuration::parameters)
        .def_ro("input", &alps::run_configuration::input)
        .def_ro("output", &alps::run_configuration::output)
        .def_ro("execution", &alps::run_configuration::execution)
        .def_ro("origins", &alps::run_configuration::origins)
        .def("save", &alps::run_configuration::save)
        .def("load", &alps::run_configuration::load);
    pyalps::mark_archive_savable(m.attr("RunConfiguration"));
    m.def("load_run_configuration", [](const std::string &filename, const std::string &schema) {
        return alps::load_run_configuration(filename, schema);
    });
    m.def("resolve_run_configuration", [](const std::string &schema, const nb::dict &parameters,
           const nb::dict &input, const nb::dict &output, const nb::dict &execution,
           const std::string &base) {
        return pyalps::resolve_run(schema, parameters, input, output, execution, base);
    });
}
