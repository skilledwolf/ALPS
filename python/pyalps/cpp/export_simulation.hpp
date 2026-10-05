// Copyright (C) 2010-2012 by Lukas Gamper
//               2026      by the ALPS collaboration
// SPDX-License-Identifier: MIT
//
// Python-owned nanobind support for downstream ALPS simulations.
// The pyalps::runtime CMake target supplies this header's include directory.
#ifndef PYALPS_EXPORT_SIMULATION_HPP
#define PYALPS_EXPORT_SIMULATION_HPP

#ifdef Py_LIMITED_API
#error "pyalps consumers must use the same per-interpreter CPython ABI as pyalps"
#endif

#include <alps/hdf5/archive.hpp>
#include <alps/mcbase.hpp>
#include "archive_savable.hpp"

#include <nanobind/nanobind.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/variant.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/shared_ptr.h>

#include <cstddef>
#include <string>
#include <utility>

NB_MAKE_OPAQUE(alps::mcbase::observable_collection_type);

namespace alps {
namespace python {

namespace nb = nanobind;

template <typename Simulation>
void export_sim_to_python(nb::module_ & module, char const * name) {
    // nanobind's type registry is shared across extension modules. Import the
    // owning pyalps modules before declaring a derived simulation so mcbase,
    // params, archive, result and observable types are already registered.
    nb::module_::import_("pyalps.ngs");
    nb::module_::import_("pyalps.hdf5");

    nb::class_<Simulation, alps::mcbase>(module, name)
        .def(nb::init<typename Simulation::parameters_type const &, std::size_t>(),
             nb::arg("parameters"), nb::arg("seed_offset") = 0)
        .def("run", [](Simulation & self, nb::object stop_callback) {
            return self.run([stop_callback = std::move(stop_callback)] {
                nb::gil_scoped_acquire gil;
                return nb::cast<bool>(stop_callback());
            });
        }, nb::arg("stop_callback"))
        .def("resultNames", &Simulation::result_names)
        .def("collectResults", [](Simulation const & self,
                                   typename Simulation::result_names_type const & names) {
            return names.empty() ? self.collect_results() : self.collect_results(names);
        }, nb::arg("names") = typename Simulation::result_names_type())
        .def("save",
             [](Simulation const & self, nb::handle archive) {
                 pyalps::with_native_archive(archive, [&](auto & native) { self.save(native); });
             })
        .def("load",
             [](Simulation & self, nb::handle archive) {
                 pyalps::with_native_archive(archive, [&](auto & native) { self.load(native); });
             });
}

}  // namespace python
}  // namespace alps

#define ALPS_EXPORT_SIM_TO_PYTHON(NAME, CLASS) \
    ::alps::python::export_sim_to_python<CLASS>(m, #NAME)

#endif
