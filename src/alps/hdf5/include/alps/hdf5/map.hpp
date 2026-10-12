/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2012 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef ALPS_NGS_HDF5_STD_MAP
#define ALPS_NGS_HDF5_STD_MAP

#include <alps/hdf5/archive.hpp>
#include <alps/ngs/cast.hpp>

#include <map>
#include <set>
#include <utility>

namespace alps {
    namespace hdf5 {

        template <typename K, typename T, typename C, typename A> void save(
              archive & ar
            , std::string const & path
            , std::map<K, T, C, A> const & value
            , std::vector<std::size_t> size = std::vector<std::size_t>()
            , std::vector<std::size_t> chunk = std::vector<std::size_t>()
            , std::vector<std::size_t> offset = std::vector<std::size_t>()
        ) {
            // Bespoke adaptation of the explicit empty-group/collection contract
            // in ALPS reconciliation 81e7ad00d. Preserve owned-group attributes.
            auto const target = ar.complete_path(path);
            if (!size.empty() || !chunk.empty() || !offset.empty()
                || target.find('@') != std::string::npos)
                throw invalid_path("maps require a complete group path" + ALPS_STACKTRACE);
            std::vector<std::string> names;
            std::set<std::string> unique;
            for (auto const & item : value) {
                auto name = ar.encode_segment(cast<std::string>(item.first));
                if (name.empty() || name == "." || name == ".."
                    || name.find('@') != std::string::npos || name.find('\0') != std::string::npos
                    || !unique.insert(name).second)
                    throw invalid_path("map keys must name distinct archive children" + ALPS_STACKTRACE);
                names.push_back(std::move(name));
            }
            ar.create_group(target);
            for (auto const & child : ar.list_children(target)) {
                auto const child_path = target + "/" + child;
                if (ar.is_group(child_path)) ar.delete_group(child_path);
                else ar.delete_data(child_path);
            }
            std::size_t index = 0;
            for (auto const & item : value)
                save(ar, target + "/" + names[index++], item.second);
        }

        template <typename K, typename T, typename C, typename A> void load(
              archive & ar
            , std::string const & path
            , std::map<K, T, C, A> & value
            , std::vector<std::size_t> chunk = std::vector<std::size_t>()
            , std::vector<std::size_t> offset = std::vector<std::size_t>()
        ) {
            if (!chunk.empty() || !offset.empty())
                throw invalid_path("maps require a complete group path" + ALPS_STACKTRACE);
            auto const target = ar.complete_path(path);
            std::map<K, T, C, A> replacement(value.key_comp(), value.get_allocator());
            for (auto const & child : ar.list_children(target)) {
                auto key = cast<K>(ar.decode_segment(child));
                auto inserted = replacement.emplace(std::move(key), T{});
                if (!inserted.second)
                    throw invalid_path("archive children decode to duplicate map keys" + ALPS_STACKTRACE);
                load(ar, target + "/" + child, inserted.first->second);
            }
            value.swap(replacement);
        }
    }
}

#endif
