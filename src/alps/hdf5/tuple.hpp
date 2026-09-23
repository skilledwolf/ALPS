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

#ifndef ALPS_NGS_HDF5_BOOST_TUPLE
#define ALPS_NGS_HDF5_BOOST_TUPLE

#include <alps/hdf5/archive.hpp>
#include <alps/ngs/detail/remove_cvr.hpp>
#include <boost/tuple/tuple.hpp>
#include <tuple>
#include <type_traits>
#include <utility>

namespace alps { namespace hdf5 {
namespace detail {

// Visit real tuple fields in index order, keeping null slots absent from the
// archive. The visitor supplies saving/loading without duplicating traversal.
template<class Types, class Tuple, class Visitor, std::size_t... I>
void visit_tuple(Tuple& value, Visitor visit, std::index_sequence<I...>) {
    auto element = [&](auto index) {
        using type = std::tuple_element_t<decltype(index)::value, Types>;
        if constexpr (!std::is_same_v<type, boost::tuples::null_type>)
            visit(index, boost::get<decltype(index)::value>(value));
    };
    (element(std::integral_constant<std::size_t, I>{}), ...);
}
} // namespace detail

template<class... T>
void save(archive& ar, std::string const& path, boost::tuple<T...> const& value,
    std::vector<std::size_t> size = {}, std::vector<std::size_t> chunk = {},
    std::vector<std::size_t> offset = {}) {
    detail::visit_tuple<std::tuple<T...>>(value, [&](auto index, auto const& element) {
        auto child = ar.complete_path(path) + "/" + std::to_string(index);
        save(ar, child, element);
        using type = typename alps::detail::remove_cvr<decltype(element)>::type;
        if (has_complex_elements<type>::value) ar.set_complex(child);
    }, std::index_sequence_for<T...>{});
}

template<class... T>
void load(archive& ar, std::string const& path, boost::tuple<T...>& value,
    std::vector<std::size_t> chunk = {}, std::vector<std::size_t> offset = {}) {
    detail::visit_tuple<std::tuple<T...>>(value, [&](auto index, auto& element) {
        load(ar, ar.complete_path(path) + "/" + std::to_string(index), element);
    }, std::index_sequence_for<T...>{});
}

}} // namespace alps::hdf5
#endif
