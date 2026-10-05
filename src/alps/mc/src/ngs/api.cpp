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

    template<class R>
    void write_results(std::map<std::string,R> const& results,
                       params const& parameters, hdf5::archive& ar, std::string const& path) {
        ar["/parameters"] << parameters;
        auto const target = ar.complete_path(path);
        ar.create_group(target);
        for (auto const& child : ar.list_children(target)) {
            auto const key = target + "/" + child;
            if (ar.is_group(key)) ar.delete_group(key);
            else ar.delete_data(key);
        }
        alea::hdf5_serializer codec(ar,target);
        for (auto const& [name,value]:results) {
            auto save=[&](auto const& result) { alea::serialize(codec,ar.encode_segment(name),result); };
            if constexpr (alea::is_alea_result<R>::value) save(value);
            else std::visit(save,value);
        }
    }
    void save_results(std::map<std::string,alea::batch_result<double>> const& results, params const& parameters,
                      hdf5::archive& ar, std::string const& path) {
        write_results(results,parameters,ar,path);
    }
    void save_results(std::map<std::string,alea::batch_result<double>> const& results, params const& parameters,
                      boost::filesystem::path const& filename, std::string const& path) {
        hdf5::archive ar(filename,"a");
        write_results(results,parameters,ar,path);
        ar.close();
    }

    void save_results(std::map<std::string,alea::result::variant_type> const& results, params const& parameters,
                      hdf5::archive& ar, std::string const& path) {
        write_results(results,parameters,ar,path);
    }
    void save_results(std::map<std::string,alea::result::variant_type> const& results, params const& parameters,
                      boost::filesystem::path const& filename, std::string const& path) {
        hdf5::archive ar(filename,"a");
        write_results(results,parameters,ar,path);
        ar.close();
    }
}
