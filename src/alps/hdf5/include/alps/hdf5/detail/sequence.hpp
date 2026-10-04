// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#pragma once

#include <alps/hdf5/archive.hpp>
#include <alps/ngs/cast.hpp>
#include <array>
#include <boost/array.hpp>

namespace alps::hdf5::detail {
// Container allocation and fixed-size checks belong here; native datatypes,
// transfers and selections belong to the HighFive-backed archive.
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
        const auto extent = hdf5::get_extent(value[0]);
        for (std::size_t i = 0; i < value.size(); ++i)
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

template<class C> void save_sequence(archive& ar, std::string const& path, C const& value,
                                    std::vector<std::size_t> size,
                                    std::vector<std::size_t> chunk,
                                    std::vector<std::size_t> offset) {
    using E = typename C::value_type;
    if constexpr (is_continuous<E>::value) {
        const auto extent = hdf5::get_extent(value);
        size.insert(size.end(), extent.begin(), extent.end());
        chunk.insert(chunk.end(), extent.begin(), extent.end());
        offset.insert(offset.end(), extent.size(), 0);
        ar.write(path, hdf5::get_pointer(value), size, chunk, offset);
    } else if (value.size() && hdf5::is_vectorizable(value)) {
        size.push_back(value.size());
        chunk.push_back(1);
        offset.push_back(0);
        for (std::size_t i = 0; i < value.size(); ++i) {
            offset.back() = i;
            save(ar, path, value[i], size, chunk, offset);
        }
    } else {
        if (path.find('@') != std::string::npos)
            throw archive_error("attributes need a native datatype");
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
        const auto children = ar.list_children(path);
        for (const auto& child : children) {
            const auto index = alps::cast<std::size_t>(child);
            if (index >= children.size() || std::to_string(index) != child)
                throw invalid_path("invalid container index: " + child);
        }
        sequence_resize<C>::apply(value, children.size());
        for (const auto& child : children)
            load(ar, ar.complete_path(path) + "/" + child,
                 value[alps::cast<std::size_t>(child)]);
        return;
    }
    const auto size = ar.extent(path);
    if (size.size() <= chunk.size()) throw archive_error("dimensions do not match");
    using E = typename C::value_type;
    if constexpr (is_continuous<E>::value) {
        hdf5::set_extent(value, {size.begin() + chunk.size(), size.end()});
        chunk.insert(chunk.end(), size.begin() + chunk.size(), size.end());
        offset.insert(offset.end(), size.size() - offset.size(), 0);
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
