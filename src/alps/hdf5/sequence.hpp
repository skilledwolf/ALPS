// Copyright (C) 2010-2011 by Lukas Gamper <gamperl@gmail.com>
// SPDX-License-Identifier: MIT
#pragma once

#include <alps/hdf5/archive.hpp>
#include <alps/ngs/cast.hpp>
#include <boost/type_traits/is_scalar.hpp>
#include <boost/type_traits/is_enum.hpp>
#include <algorithm>
#include <iterator>
#include <utility>

namespace alps { namespace hdf5 { namespace detail {
// Indexed sequences share their archive representation. FixedExtent selects
// array shape validation without resizing; packed vector<bool> keeps its adapter.
template<class Sequence, bool FixedExtent = false> struct sequence_get_extent {
    static std::vector<std::size_t> apply(Sequence const& value) {
        using alps::hdf5::get_extent;
        std::vector<std::size_t> result(1, value.size());
        if (!value.size()) return result;
        auto first = get_extent(value[0]);
        if (!FixedExtent && !boost::is_scalar<typename Sequence::value_type>::value)
            for (std::size_t i = 1; i < value.size(); ++i)
                if (get_extent(value[i]) != first)
                    throw archive_error("no rectangular matrix" + ALPS_STACKTRACE);
        result.insert(result.end(), first.begin(), first.end());
        return result;
    }
};

template<class Sequence, bool FixedExtent = false> struct sequence_set_extent {
    static void apply(Sequence& value, std::vector<std::size_t> const& extent) {
        using alps::hdf5::set_extent;
        using element_type = typename Sequence::value_type;
        using scalar = typename scalar_type<element_type>::type;
        if constexpr (!FixedExtent) {
            if (extent.empty()) throw archive_error("dimensions do not match" + ALPS_STACKTRACE);
            value.resize(extent[0]);
        }
        if (extent.size() > 1) {
            std::vector<std::size_t> inner(extent.begin() + 1, extent.end());
            for (std::size_t i = 0; i < value.size(); ++i) set_extent(value[i], inner);
        } else if constexpr (FixedExtent) {
            if (extent.empty() && !boost::is_same<scalar, element_type>::value)
                throw archive_error("dimensions do not match" + ALPS_STACKTRACE);
        } else if ((!boost::is_enum<element_type>::value && !boost::is_same<scalar, element_type>::value) ||
                   (boost::is_enum<element_type>::value && is_continuous<element_type>::value && sizeof(element_type) != sizeof(scalar))) {
            throw archive_error("dimensions do not match" + ALPS_STACKTRACE);
        }
    }
};

template<class Sequence, bool FixedExtent = false> struct sequence_is_vectorizable {
    static bool apply(Sequence const& value) {
        using alps::hdf5::get_extent;
        using alps::hdf5::is_vectorizable;
        if (!value.size() || (FixedExtent && is_continuous<Sequence>::value)) return true;
        if (!is_vectorizable(value[0])) return false;
        auto first = get_extent(value[0]);
        if (!boost::is_scalar<typename Sequence::value_type>::value)
            for (std::size_t i = 0; i < value.size(); ++i)
                if (!is_vectorizable(value[i]) || get_extent(value[i]) != first)
                    return false;
        return true;
    }
};

template<class Sequence> struct sequence_get_pointer {
    static auto apply(Sequence& value) {
        using alps::hdf5::get_pointer;
        return get_pointer(value[0]);
    }
};

template<class Sequence> void save_sequence(archive& ar, std::string const& path, Sequence const& value,
    std::vector<std::size_t> size, std::vector<std::size_t> chunk, std::vector<std::size_t> offset) {
    using alps::cast;
    using alps::hdf5::get_extent;
    using alps::hdf5::set_extent;
    using alps::hdf5::get_pointer;
    using alps::hdf5::is_vectorizable;
    if (ar.is_group(path))
        ar.delete_group(path);
    if (is_continuous<typename Sequence::value_type>::value && value.size() == 0)
        ar.write(path, static_cast<typename scalar_type<Sequence >::type const *>(NULL), std::vector<std::size_t>());
    else if (is_continuous<typename Sequence::value_type>::value) {
        std::vector<std::size_t> extent(get_extent(value));
        std::copy(extent.begin(), extent.end(), std::back_inserter(size));
        std::copy(extent.begin(), extent.end(), std::back_inserter(chunk));
        std::fill_n(std::back_inserter(offset), extent.size(), 0);
        ar.write(path, get_pointer(value), size, chunk, offset);
    } else if (value.size() == 0)
        ar.write(path, static_cast<int const *>(NULL), std::vector<std::size_t>());
    else if (is_vectorizable(value)) {
        size.push_back(value.size());
        chunk.push_back(1);
        offset.push_back(0);
        for (std::size_t i = 0; i < value.size(); ++i) {
            offset.back() = i;
            save(ar, path, value[i], size, chunk, offset);
        }
    } else {
        if (path.find_last_of('@') == std::string::npos && ar.is_data(path))
            ar.delete_data(path);
        else if (path.find_last_of('@') != std::string::npos && ar.is_attribute(path))
            ar.delete_attribute(path);
        for (std::size_t i = 0; i < value.size(); ++i)
            save(ar, ar.complete_path(path) + "/" + cast<std::string>(i), value[i]);
    }
}

template<class Sequence, bool FixedExtent = false> void load_sequence(archive& ar, std::string const& path, Sequence& value,
    std::vector<std::size_t> chunk, std::vector<std::size_t> offset) {
    using alps::cast;
    using alps::hdf5::get_extent;
    using alps::hdf5::set_extent;
    using alps::hdf5::get_pointer;
    using alps::hdf5::is_vectorizable;
    if (ar.is_group(path)) {
        std::vector<std::string> children = ar.list_children(path);
        if constexpr (FixedExtent) {
            if (children.size() != value.size())
                throw invalid_path("size does not match: " + path + ALPS_STACKTRACE);
        } else {
            value.resize(children.size());
        }
        for (typename std::vector<std::string>::const_iterator it = children.begin(); it != children.end(); ++it)
           load(ar, ar.complete_path(path) + "/" + *it, value[cast<std::size_t>(*it)]);
    } else {
        if (ar.is_complex(path) != has_complex_elements<typename Sequence::value_type>::value)
            throw archive_error("no complex value in archive" + ALPS_STACKTRACE);
        std::vector<std::size_t> size(ar.extent(path));
        if constexpr (FixedExtent) {
            if (!size.empty() && (chunk.size() >= size.size() ||
                (value.size() != size[chunk.size()] &&
                 (is_continuous<typename Sequence::value_type>::value || size[chunk.size()] > 0))))
                throw archive_error("dimensions do not match" + ALPS_STACKTRACE);
        } else {
            if (size.empty()) throw archive_error("invalid dimensions" + ALPS_STACKTRACE);
            if (size[0] == 0) {
                value.resize(0);
                return;
            }
        }
        if (is_continuous<typename Sequence::value_type>::value) {
            set_extent(value, std::vector<std::size_t>(size.begin() + chunk.size(), size.end()));
            if (value.size()) {
                std::copy(size.begin() + chunk.size(), size.end(), std::back_inserter(chunk));
                std::fill_n(std::back_inserter(offset), size.size() - offset.size(), 0);
                ar.read(path, get_pointer(value), chunk, offset);
            }
        } else {
            if (chunk.size() >= size.size())
                throw archive_error("dimensions do not match" + ALPS_STACKTRACE);
            if constexpr (FixedExtent)
                set_extent(value, std::vector<std::size_t>(1, size[chunk.size()]));
            else
                value.resize(size[chunk.size()]);
            chunk.push_back(1);
            offset.push_back(0);
            for (std::size_t i = 0; i < value.size(); ++i) {
                offset.back() = i;
                load(ar, path, value[i], chunk, offset);
            }
        }
    }
}

}}} // namespace alps::hdf5::detail
