// Copyright (C) 2010–2012 Lukas Gamper; 2026 ALPS Collaboration.
// SPDX-License-Identifier: MIT
#pragma once
#include <alps/hdf5/detail/sequence.hpp>
#include <array>

namespace alps::hdf5 {
template<typename T, std::size_t N> struct scalar_type<std::array<T, N>> : scalar_type<T> {};
template<typename T, std::size_t N> struct is_continuous<std::array<T, N>> : is_continuous<T> {};
template<typename T, std::size_t N> struct is_continuous<std::array<T, N> const> : is_continuous<T> {};
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
