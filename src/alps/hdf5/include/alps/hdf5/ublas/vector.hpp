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

// Adapted from ALPS reconciliation 2cc35a558, corresponding sequence adapter.
// Keep the current C backend, legacy empty schema and complex type markers.
#pragma once
#include <alps/hdf5/detail/sequence.hpp>
#include <boost/numeric/ublas/vector.hpp>

namespace alps::hdf5 {
template<typename T, typename A> struct scalar_type<boost::numeric::ublas::vector<T, A>> : scalar_type<T> {};
template<typename T, typename A> struct has_complex_elements<boost::numeric::ublas::vector<T, A>> : has_complex_elements<T> {};
namespace detail {
template<typename T, typename A> struct get_extent<boost::numeric::ublas::vector<T, A>> : sequence_extent<boost::numeric::ublas::vector<T, A>> {};
template<typename T, typename A> struct set_extent<boost::numeric::ublas::vector<T, A>> : sequence_set_extent<boost::numeric::ublas::vector<T, A>> {};
template<typename T, typename A> struct is_vectorizable<boost::numeric::ublas::vector<T, A>> : sequence_vectorizable<boost::numeric::ublas::vector<T, A>> {};
template<typename T, typename A> struct get_pointer<boost::numeric::ublas::vector<T, A>> : sequence_pointer<boost::numeric::ublas::vector<T, A>> {};
template<typename T, typename A> struct get_pointer<boost::numeric::ublas::vector<T, A> const> : sequence_pointer<boost::numeric::ublas::vector<T, A> const> {};
} // namespace detail

template<typename T, typename A> void save(archive& ar, std::string const& path, boost::numeric::ublas::vector<T, A> const& value,
    std::vector<std::size_t> size = {}, std::vector<std::size_t> chunk = {},
    std::vector<std::size_t> offset = {}) {
    detail::save_sequence(ar, path, value, std::move(size), std::move(chunk), std::move(offset));
}
template<typename T, typename A> void load(archive& ar, std::string const& path, boost::numeric::ublas::vector<T, A>& value,
    std::vector<std::size_t> chunk = {}, std::vector<std::size_t> offset = {}) {
    detail::load_sequence(ar, path, value, std::move(chunk), std::move(offset));
}
} // namespace alps::hdf5
