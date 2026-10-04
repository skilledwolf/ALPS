// Copyright (C) 2010–2012 Lukas Gamper; 2026 ALPS Collaboration.
// SPDX-License-Identifier: MIT
#pragma once
#include <alps/hdf5/detail/sequence.hpp>
#include <vector>

namespace alps::hdf5 {
template<typename T, typename A> struct scalar_type<std::vector<T, A>> : scalar_type<T> {};
template<typename T, typename A> struct is_content_continuous<std::vector<T, A>> : is_continuous<T> {};
template<typename A> struct is_content_continuous<std::vector<bool, A>> : boost::false_type {};
namespace detail {
template<typename T, typename A> struct get_extent<std::vector<T, A>> : sequence_extent<std::vector<T, A>> {};
template<typename T, typename A> struct set_extent<std::vector<T, A>> : sequence_set_extent<std::vector<T, A>> {};
template<typename T, typename A> struct is_vectorizable<std::vector<T, A>> : sequence_vectorizable<std::vector<T, A>> {};
template<typename T, typename A> struct get_pointer<std::vector<T, A>> : sequence_pointer<std::vector<T, A>> {};
template<typename T, typename A> struct get_pointer<std::vector<T, A> const> : sequence_pointer<std::vector<T, A> const> {};
template<typename A> struct is_vectorizable<std::vector<bool, A>> {
    static bool apply(std::vector<bool, A> const&) { return true; }
};
template<typename A> struct set_extent<std::vector<bool, A>> {
    static void apply(std::vector<bool, A>& value, std::vector<std::size_t> const& extent) {
        if (extent.size() != 1) throw archive_error("dimensions do not match");
        value.resize(extent[0]);
    }
};
template<typename A> struct get_pointer<std::vector<bool, A>> {
    static bool* apply(std::vector<bool, A>&) {
        throw archive_error("packed Boolean vectors have no native pointer");
    }
};
template<typename A> struct get_pointer<std::vector<bool, A> const> {
    static bool const* apply(std::vector<bool, A> const&) {
        throw archive_error("packed Boolean vectors have no native pointer");
    }
};
} // namespace detail

template<typename T, typename A> void save(archive& ar, std::string const& path, std::vector<T, A> const& value,
    std::vector<std::size_t> size = {}, std::vector<std::size_t> chunk = {},
    std::vector<std::size_t> offset = {}) {
    detail::save_sequence(ar, path, value, std::move(size), std::move(chunk), std::move(offset));
}
template<typename T, typename A> void load(archive& ar, std::string const& path, std::vector<T, A>& value,
    std::vector<std::size_t> chunk = {}, std::vector<std::size_t> offset = {}) {
    detail::load_sequence(ar, path, value, std::move(chunk), std::move(offset));
}

// std::vector<bool> packs bits; unpack once for one native transfer rather than
// making a separate HDF5 selection for every element.
template<typename A> void save(archive& ar, std::string const& path, std::vector<bool, A> const& value,
    std::vector<std::size_t> size = {}, std::vector<std::size_t> chunk = {},
    std::vector<std::size_t> offset = {}) {
    auto buffer = std::make_unique<bool[]>(value.size());
    for (std::size_t i = 0; i < value.size(); ++i) buffer[i] = value[i];
    size.push_back(value.size());
    chunk.push_back(value.size());
    offset.push_back(0);
    ar.write(path, buffer.get(), size, chunk, offset);
}
template<typename A> void load(archive& ar, std::string const& path, std::vector<bool, A>& value,
    std::vector<std::size_t> chunk = {}, std::vector<std::size_t> offset = {}) {
    const auto size = ar.extent(path);
    if (ar.is_group(path) || size.size() != chunk.size() + 1)
        throw archive_error("dimensions do not match");
    const auto count = size.back();
    auto buffer = std::make_unique<bool[]>(count);
    chunk.push_back(count);
    offset.push_back(0);
    ar.read(path, buffer.get(), chunk, offset);
    value.assign(buffer.get(), buffer.get() + count);
}
} // namespace alps::hdf5
