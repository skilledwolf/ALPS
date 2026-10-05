// Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>
//               Matthias Troyer <troyer@comp-phys.org>
//               2026       by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
#include <nanobind/nanobind.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/string.h>
#include <alps/ngs/api.hpp>
#include <alps/params.hpp>
#include <alps/hdf5/archive.hpp>
#include <alps/mcbase.hpp>
#include "archive_savable.hpp"
namespace nb = nanobind;
namespace alps {
    namespace detail {
        void save_results_export(mcbase::results_type const & res, params const & par, nb::handle object, std::string const & path) {
            pyalps::with_native_archive(object, [&](auto & ar) {
                alps::save_results(res, par, ar, path);
            });
        }
    }
}
NB_MODULE(pyngsapi_c, m) {
    m.def("collectResults", [](alps::mcbase const & sim) {
        return sim.collect_results();
    });
    m.def("saveResults", &alps::detail::save_results_export);
}
