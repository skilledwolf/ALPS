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
#include <alps/alea/hdf5.hpp>

namespace alps {

    void save_results(std::map<std::string, alps::alea::batch_result<double>> const & results,
                      params const & params, boost::filesystem::path const & filename,
                      std::string const & path) {
        hdf5::archive ar(filename, "a");
        save_results(results, params, ar, path);
        ar.close();
    }

    void save_results(std::map<std::string, alps::alea::batch_result<double>> const & results,
                      params const & params, hdf5::archive & ar, std::string const & path) {
        ar["/parameters"] << params;
        auto const target = ar.complete_path(path);
        ar.create_group(target);
        for (auto const& child : ar.list_children(target)) {
            auto const key = target + "/" + child;
            if (ar.is_group(key)) ar.delete_group(key);
            else ar.delete_data(key);
        }
        alps::alea::hdf5_serializer serializer(ar, target);
        for (auto const& entry : results)
            alps::alea::serialize(serializer, ar.encode_segment(entry.first), entry.second);
    }

}
