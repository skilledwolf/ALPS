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

// Adapted from ALPS reconciliation 2cc35a558, corresponding sequence adapter.
// Keep the current C backend, legacy empty schema and complex type markers.
#pragma once
#include <alps/hdf5/detail/sequence.hpp>
#include <array>

namespace alps::hdf5 {
template<typename T, std::size_t N> struct scalar_type<std::array<T, N>> : scalar_type<T> {};
template<typename T, std::size_t N> struct is_continuous<std::array<T, N>> : is_continuous<T> {};
template<typename T, std::size_t N> struct is_continuous<std::array<T, N> const> : is_continuous<T> {};
template<typename T, std::size_t N> struct has_complex_elements<std::array<T, N>> : has_complex_elements<T> {};
namespace detail {
template<typename T, std::size_t N> struct get_extent<std::array<T, N>> : sequence_extent<std::array<T, N>> {};
template<typename T, std::size_t N> struct set_extent<std::array<T, N>> : sequence_set_extent<std::array<T, N>> {};
template<typename T, std::size_t N> struct is_vectorizable<std::array<T, N>> : sequence_vectorizable<std::array<T, N>> {};
template<typename T, std::size_t N> struct get_pointer<std::array<T, N>> : sequence_pointer<std::array<T, N>> {};
template<typename T, std::size_t N> struct get_pointer<std::array<T, N> const> : sequence_pointer<std::array<T, N> const> {};
} // namespace detail

template<typename T, std::size_t N> void save(archive& ar, std::string const& path, std::array<T, N> const& value,
    std::vector<std::size_t> size = {}, std::vector<std::size_t> chunk = {},
    std::vector<std::size_t> offset = {}) {
    detail::save_sequence(ar, path, value, std::move(size), std::move(chunk), std::move(offset));
}
template<typename T, std::size_t N> void load(archive& ar, std::string const& path, std::array<T, N>& value,
    std::vector<std::size_t> chunk = {}, std::vector<std::size_t> offset = {}) {
    detail::load_sequence(ar, path, value, std::move(chunk), std::move(offset));
}
} // namespace alps::hdf5
