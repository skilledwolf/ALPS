// Copyright (C) 2010 - 2013 by Lukas Gamper <gamperl@gmail.com>
//               Matthias Troyer <troyer@comp-phys.org>
//               2026       by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <alps/ngs/random01.hpp>
#include "archive_savable.hpp"
#include <alps/hdf5/archive.hpp>
namespace nb = nanobind;
NB_MODULE(pyngsrandom01_c, m) {
    nb::class_<alps::random01>(m, "random01")
        .def(nb::init<int, std::string const&>(), nb::arg("seed") = 42, nb::arg("name") = "mt19937")
        .def_prop_ro("name", &alps::random01::name)
        .def("__deepcopy__",
             // copy.deepcopy() passes (self, memo); memo is unused.
             [](alps::random01 const & self, nb::handle /*memo*/) {
                 return alps::random01(self);
             })
        .def("__call__",
             static_cast<alps::random01::result_type (alps::random01::*)()>(
                 &alps::random01::operator()))
        .def("save", [](alps::random01 const & self, nb::handle ar) { pyalps::with_native_archive(ar, [&](auto & native) { self.save(native); }); })
        .def("load", [](alps::random01 & self, nb::handle ar) { pyalps::with_native_archive(ar, [&](auto & native) { self.load(native); }); });
    pyalps::mark_archive_savable(m.attr("random01"));
}
