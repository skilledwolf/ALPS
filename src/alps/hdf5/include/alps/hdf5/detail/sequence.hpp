// Copyright (C) 2010–2012 Lukas Gamper; 2026 ALPS Collaboration.
// SPDX-License-Identifier: MIT
// Adapted from ALPS reconciliation 2cc35a558, detail/sequence.hpp and the
// sequence adapters it consolidates. Preserve the C backend's legacy H5S_NULL
// empty values and complex markers; validate selections before indexing.
#pragma once

#include <alps/hdf5/archive.hpp>
#include <algorithm>
#include <array>
#include <boost/array.hpp>
#include <memory>
#include <set>
#include <utility>

namespace alps::hdf5::detail {

template<class C> struct sequence_resize {
    static void apply(C& value, std::size_t size) { value.resize(size); }
};
template<class T, std::size_t N> struct sequence_resize<std::array<T, N>> {
    static void apply(std::array<T, N>&, std::size_t size) {
        if (size != N) throw archive_error("dimensions do not match");
    }
};
template<class T, std::size_t N> struct sequence_resize<boost::array<T, N>> {
    static void apply(boost::array<T, N>&, std::size_t size) {
        if (size != N) throw archive_error("dimensions do not match");
    }
};

template<class C> struct sequence_extent {
    static std::vector<std::size_t> apply(C const& value) {
        using E = typename C::value_type;
        std::vector<std::size_t> result{value.size()}, element;
        if (value.size()) {
            element = hdf5::get_extent(value[0]);
            for (std::size_t i = 1; i < value.size(); ++i)
                if (element != hdf5::get_extent(value[i]))
                    throw archive_error("container is not rectangular");
        } else if constexpr (is_continuous<E>::value) {
            element = hdf5::get_extent(E{});
        }
        result.insert(result.end(), element.begin(), element.end());
        return result;
    }
};

template<class C> struct sequence_set_extent {
    static void apply(C& value, std::vector<std::size_t> const& extent) {
        if (extent.empty()) throw archive_error("dimensions do not match");
        using E = typename C::value_type;
        const std::vector<std::size_t> element(extent.begin() + 1, extent.end());
        // Validate fixed inner dimensions even when the outer sequence is empty.
        if constexpr (is_continuous<E>::value) {
            if (element != hdf5::get_extent(E{}))
                throw archive_error("dimensions do not match");
        }
        sequence_resize<C>::apply(value, extent[0]);
        for (std::size_t i = 0; i < value.size(); ++i)
            hdf5::set_extent(value[i], element);
    }
};

template<class C> struct sequence_vectorizable {
    static bool apply(C const& value) {
        if (!value.size()) return true;
        if (!hdf5::is_vectorizable(value[0])) return false;
        const auto extent = hdf5::get_extent(value[0]);
        for (std::size_t i = 1; i < value.size(); ++i)
            if (!hdf5::is_vectorizable(value[i]) || hdf5::get_extent(value[i]) != extent)
                return false;
        return true;
    }
};

template<class C> struct sequence_pointer {
    static typename scalar_type<C>::type* apply(C& value) {
        return value.size() ? hdf5::get_pointer(value[0]) : nullptr;
    }
};
template<class C> struct sequence_pointer<C const> {
    static typename scalar_type<C>::type const* apply(C const& value) {
        return value.size() ? hdf5::get_pointer(value[0]) : nullptr;
    }
};

inline void validate_sequence_write_prefix(std::vector<std::size_t> const& size,
                                           std::vector<std::size_t> const& chunk,
                                           std::vector<std::size_t> const& offset) {
    if (size.size() != chunk.size() || size.size() != offset.size())
        throw archive_error("dimensions do not match");
    for (std::size_t i = 0; i < size.size(); ++i)
        // A sequence owns one element of each enclosing container, not a
        // buffer spanning multiple outer elements.
        if (chunk[i] != 1 || offset[i] >= size[i])
            throw archive_error("invalid sequence selection");
}

inline void validate_sequence_read_prefix(std::vector<std::size_t> const& size,
                                          std::vector<std::size_t> const& chunk,
                                          std::vector<std::size_t> const& offset) {
    if (chunk.size() != offset.size() || size.size() <= chunk.size())
        throw archive_error("dimensions do not match");
    for (std::size_t i = 0; i < chunk.size(); ++i)
        if (chunk[i] != 1 || offset[i] >= size[i])
            throw archive_error("invalid sequence selection");
}

inline void validate_sequence_children(std::vector<std::string> const& children) {
    // Compare names with the complete expected index set. Parsing arbitrary
    // archive strings as integers can accept signs/leading zeroes or overflow.
    std::set<std::string> names(children.begin(), children.end());
    for (std::size_t i = 0; i < children.size(); ++i)
        if (names.count(std::to_string(i)) == 0)
            throw invalid_path("container indices must be dense and zero based");
}

template<class C> void save_sequence(archive& ar, std::string const& path, C const& value,
                                    std::vector<std::size_t> size,
                                    std::vector<std::size_t> chunk,
                                    std::vector<std::size_t> offset) {
    using E = typename C::value_type;
    validate_sequence_write_prefix(size, chunk, offset);
    if (value.size() == 0 && size.empty()) {
        // Keep the historical top-level null datatype, including int for
        // sequences of non-contiguous values. Nested empties retain their rank.
        if constexpr (is_continuous<E>::value)
            ar.write(path, static_cast<typename scalar_type<C>::type const*>(nullptr), {});
        else ar.write(path, static_cast<int const*>(nullptr), {});
    } else if constexpr (is_continuous<E>::value) {
        const auto extent = hdf5::get_extent(value);
        size.insert(size.end(), extent.begin(), extent.end());
        chunk.insert(chunk.end(), extent.begin(), extent.end());
        offset.insert(offset.end(), extent.size(), 0);
        ar.write(path, hdf5::get_pointer(value), size, chunk, offset);
    } else if (value.size() == 0) {
        size.push_back(0);
        chunk.push_back(0);
        offset.push_back(0);
        ar.write(path, static_cast<int const*>(nullptr), size, chunk, offset);
    } else if (hdf5::is_vectorizable(value)) {
        size.push_back(value.size());
        chunk.push_back(1);
        offset.push_back(0);
        for (std::size_t i = 0; i < value.size(); ++i) {
            offset.back() = i;
            save(ar, path, value[i], size, chunk, offset);
        }
    } else {
        if (!size.empty() || path.find('@') != std::string::npos)
            throw archive_error("ragged sequences require an unselected group path");
        if (ar.is_group(path)) ar.delete_group(path);
        if (ar.is_data(path)) ar.delete_data(path);
        ar.create_group(path);
        for (std::size_t i = 0; i < value.size(); ++i)
            save(ar, ar.complete_path(path) + "/" + std::to_string(i), value[i]);
    }
}

template<class C> void load_sequence(archive& ar, std::string const& path, C& value,
                                    std::vector<std::size_t> chunk,
                                    std::vector<std::size_t> offset) {
    if (ar.is_group(path)) {
        if (!chunk.empty() || !offset.empty())
            throw archive_error("ragged sequences require an unselected group path");
        const auto children = ar.list_children(path);
        validate_sequence_children(children);
        sequence_resize<C>::apply(value, children.size());
        for (std::size_t i = 0; i < children.size(); ++i)
            load(ar, ar.complete_path(path) + "/" + std::to_string(i), value[i]);
        return;
    }
    using E = typename C::value_type;
    if (ar.is_complex(path) != has_complex_elements<E>::value)
        throw archive_error("no complex value in archive");
    if (ar.is_scalar(path)) throw archive_error("dimensions do not match");
    const auto size = ar.extent(path);
    validate_sequence_read_prefix(size, chunk, offset);
    if (ar.is_null(path)) {
        sequence_resize<C>::apply(value, 0);
        return;
    }
    if constexpr (is_continuous<E>::value) {
        const std::vector<std::size_t> extent(size.begin() + chunk.size(), size.end());
        hdf5::set_extent(value, extent);
        if (std::find(extent.begin(), extent.end(), 0) != extent.end()) return;
        chunk.insert(chunk.end(), extent.begin(), extent.end());
        offset.insert(offset.end(), extent.size(), 0);
        ar.read(path, hdf5::get_pointer(value), chunk, offset);
    } else {
        sequence_resize<C>::apply(value, size[chunk.size()]);
        chunk.push_back(1);
        offset.push_back(0);
        for (std::size_t i = 0; i < value.size(); ++i) {
            offset.back() = i;
            load(ar, path, value[i], chunk, offset);
        }
    }
}
} // namespace alps::hdf5::detail
