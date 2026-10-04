// Copyright (C) 2010–2012 Lukas Gamper; 2026 ALPS Collaboration.
// SPDX-License-Identifier: MIT
#pragma once
#include <alps/hdf5/detail/sequence.hpp>
#include <valarray>

namespace alps::hdf5 {
template<typename T> struct scalar_type<std::valarray<T>> : scalar_type<T> {};
template<typename T> struct is_content_continuous<std::valarray<T>> : is_continuous<T> {};
namespace detail {
template<typename T> struct get_extent<std::valarray<T>> : sequence_extent<std::valarray<T>> {};
template<typename T> struct set_extent<std::valarray<T>> : sequence_set_extent<std::valarray<T>> {};
template<typename T> struct is_vectorizable<std::valarray<T>> : sequence_vectorizable<std::valarray<T>> {};
template<typename T> struct get_pointer<std::valarray<T>> : sequence_pointer<std::valarray<T>> {};
template<typename T> struct get_pointer<std::valarray<T> const> : sequence_pointer<std::valarray<T> const> {};
} // namespace detail

template<typename T> void save(archive& ar, std::string const& path, std::valarray<T> const& value,
    std::vector<std::size_t> size = {}, std::vector<std::size_t> chunk = {},
    std::vector<std::size_t> offset = {}) {
    detail::save_sequence(ar, path, value, std::move(size), std::move(chunk), std::move(offset));
}
template<typename T> void load(archive& ar, std::string const& path, std::valarray<T>& value,
    std::vector<std::size_t> chunk = {}, std::vector<std::size_t> offset = {}) {
    detail::load_sequence(ar, path, value, std::move(chunk), std::move(offset));
}
} // namespace alps::hdf5
