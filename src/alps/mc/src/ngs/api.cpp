/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <alps/ngs/api.hpp>
#include <alps/hdf5/archive.hpp>

#include <boost/filesystem.hpp>

namespace alps {

    namespace detail {
        template<typename R, typename P> void save_results_impl(R const & results, P const & params, boost::filesystem::path const & filename, std::string const & path) {
            // Adapted from ALPS reconciliation 81e7ad00d, ngs/api.cpp.
            // An empty result set replaces its group too; unrelated groups stay.
            hdf5::archive ar(filename.string(), "a");
            auto const target = ar.complete_path(path);
            if (target == "/" || target == "/parameters" || target.compare(0, 12, "/parameters/") == 0)
                throw hdf5::invalid_path("result path overlaps the parameters group: " + target);
            ar["/parameters"] << params;
            ar[target] << results;
            ar.close();
        }
    }

    #ifdef ALPS_NGS_USE_NEW_ALEA

        void save_results(alps::accumulator::result_set const & results, params const & params, boost::filesystem::path const & filename, std::string const & path) {
            detail::save_results_impl(results, params, filename, path);
        }

        void save_results(alps::accumulator::accumulator_set const & observables, params const & params, boost::filesystem::path const & filename, std::string const & path) {
            detail::save_results_impl(observables, params, filename, path);
        }

    #endif

    void save_results(mcresults const & results, params const & params, boost::filesystem::path const & filename, std::string const & path) {
        detail::save_results_impl(results, params, filename, path);
    }

    void save_results(mcobservables const & observables, params const & params, boost::filesystem::path const & filename, std::string const & path) {
        detail::save_results_impl(observables, params, filename, path);
    }

}
