// Copyright (C) 2010–2012 Lukas Gamper; 2026 ALPS Collaboration.
// SPDX-License-Identifier: MIT
#pragma once
#include <alps/hdf5/detail/sequence.hpp>
#include <boost/array.hpp>

namespace alps::hdf5 {
template<typename T, std::size_t N> struct scalar_type<boost::array<T, N>> : scalar_type<T> {};
template<typename T, std::size_t N> struct is_continuous<boost::array<T, N>> : is_continuous<T> {};
template<typename T, std::size_t N> struct is_continuous<boost::array<T, N> const> : is_continuous<T> {};
namespace detail {
template<typename T, std::size_t N> struct get_extent<boost::array<T, N>> : sequence_extent<boost::array<T, N>> {};
template<typename T, std::size_t N> struct set_extent<boost::array<T, N>> : sequence_set_extent<boost::array<T, N>> {};
template<typename T, std::size_t N> struct is_vectorizable<boost::array<T, N>> : sequence_vectorizable<boost::array<T, N>> {};
template<typename T, std::size_t N> struct get_pointer<boost::array<T, N>> : sequence_pointer<boost::array<T, N>> {};
template<typename T, std::size_t N> struct get_pointer<boost::array<T, N> const> : sequence_pointer<boost::array<T, N> const> {};
} // namespace detail

template<typename T, std::size_t N> void save(archive& ar, std::string const& path, boost::array<T, N> const& value,
    std::vector<std::size_t> size = {}, std::vector<std::size_t> chunk = {},
    std::vector<std::size_t> offset = {}) {
    detail::save_sequence(ar, path, value, std::move(size), std::move(chunk), std::move(offset));
}
template<typename T, std::size_t N> void load(archive& ar, std::string const& path, boost::array<T, N>& value,
    std::vector<std::size_t> chunk = {}, std::vector<std::size_t> offset = {}) {
    detail::load_sequence(ar, path, value, std::move(chunk), std::move(offset));
}
} // namespace alps::hdf5
