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

#ifndef ALPS_NGS_HDF5_STD_PAIR
#define ALPS_NGS_HDF5_STD_PAIR

#include <alps/hdf5/archive.hpp>
#include <alps/ngs/cast.hpp>
#include <alps/ngs/detail/remove_cvr.hpp>

#include <utility>
#include <algorithm>
#include <numeric>

namespace alps {
    namespace hdf5 {

        template <typename T, typename U> void save(
              archive & ar
            , std::string const & path
            , std::pair<T, U> const & value
            , std::vector<std::size_t> size = std::vector<std::size_t>()
            , std::vector<std::size_t> chunk = std::vector<std::size_t>()
            , std::vector<std::size_t> offset = std::vector<std::size_t>()
        ) {
            save(ar, ar.complete_path(path) + "/0", value.first);
            if (has_complex_elements<typename alps::detail::remove_cvr<T>::type>::value)
                ar.set_complex(ar.complete_path(path) + "/0");
            save(ar, ar.complete_path(path) + "/1", value.second);
            if (has_complex_elements<typename alps::detail::remove_cvr<U>::type>::value)
                ar.set_complex(ar.complete_path(path) + "/1");
        }

        template <typename T, typename U> void load(
              archive & ar
            , std::string const & path
            , std::pair<T, U> & value
            , std::vector<std::size_t> chunk = std::vector<std::size_t>()
            , std::vector<std::size_t> offset = std::vector<std::size_t>()
        ) {
            try {
                load(ar, ar.complete_path(path) + "/0", value.first);
                load(ar, ar.complete_path(path) + "/1", value.second);
            } catch (path_not_found exc) {
                load(ar, ar.complete_path(path) + "/first", value.first);
                load(ar, ar.complete_path(path) + "/second", value.second);
            }
        }

        template<typename T> struct scalar_type<std::pair<T*, std::vector<std::size_t>>> {
            using type = typename scalar_type<typename alps::detail::remove_cvr<T>::type>::type;
        };
        template<typename T> struct is_content_continuous<std::pair<T*, std::vector<std::size_t>>>
            : is_continuous<T> {};
        template<typename T> struct has_complex_elements<std::pair<T*, std::vector<std::size_t>>>
            : has_complex_elements<typename alps::detail::remove_cvr<T>::type> {};

        namespace detail {
            inline std::size_t array_size(std::vector<std::size_t> const& shape) {
                return shape.empty() ? 0 : std::accumulate(shape.begin(), shape.end(), std::size_t(1), std::multiplies<std::size_t>());
            }

            // Traverse a flat buffer in the archive's row-major coordinate
            // order. Both dataset hyperslabs and ragged groups use these indices.
            template<class Visitor>
            void visit_array(std::vector<std::size_t> const& shape,
                std::vector<std::size_t> offset, Visitor visit) {
                auto prefix = offset.size();
                offset.resize(prefix + shape.size(), 0);
                auto count = array_size(shape);
                for (std::size_t i = 0; i < count; ++i) {
                    visit(i, offset);
                    for (std::size_t d = shape.size(); d-- > 0;) {
                        if (++offset[prefix + d] < shape[d]) break;
                        offset[prefix + d] = 0;
                    }
                }
            }

            inline std::string array_path(std::string path, std::vector<std::size_t> const& index) {
                for (auto i : index) path += "/" + std::to_string(i);
                return path;
            }

            template<typename T> struct get_extent<std::pair<T*, std::vector<std::size_t>>> {
                static std::vector<std::size_t> apply(std::pair<T*, std::vector<std::size_t>> const& value) {
                    using alps::hdf5::get_extent;
                    auto extent = value.second;
                    auto count = array_size(extent);
                    if (!count) return extent;
                    auto inner = get_extent(*value.first);
                    if (!is_continuous<T>::value)
                        for (std::size_t i = 1; i < count; ++i)
                            if (get_extent(value.first[i]) != inner)
                                throw archive_error("no rectangular matrix" + ALPS_STACKTRACE);
                    extent.insert(extent.end(), inner.begin(), inner.end());
                    return extent;
                }
            };

            template<typename T> struct set_extent<std::pair<T*, std::vector<std::size_t>>> {
                static void apply(std::pair<T*, std::vector<std::size_t>>& value, std::vector<std::size_t> const& size) {
                    using alps::hdf5::set_extent;
                    if (value.second.size() > size.size() || !std::equal(value.second.begin(), value.second.end(), size.begin()))
                        throw archive_error("invalid data size" + ALPS_STACKTRACE);
                    if (!is_continuous<T>::value && value.second.size() < size.size()) {
                        std::vector<std::size_t> inner(size.begin() + value.second.size(), size.end());
                        for (std::size_t i = 0; i < array_size(value.second); ++i) set_extent(value.first[i], inner);
                    }
                }
            };

            template<typename T> struct is_vectorizable<std::pair<T*, std::vector<std::size_t>>> {
                static bool apply(std::pair<T*, std::vector<std::size_t>> const& value) {
                    using alps::hdf5::get_extent;
                    using alps::hdf5::is_vectorizable;
                    auto count = array_size(value.second);
                    if (!count) return true;
                    if (!is_vectorizable(value.first[0])) return false;
                    auto inner = get_extent(value.first[0]);
                    for (std::size_t i = 1; i < count; ++i)
                        if (!is_vectorizable(value.first[i]) || get_extent(value.first[i]) != inner) return false;
                    return true;
                }
            };

            template<typename T> struct get_pointer<std::pair<T*, std::vector<std::size_t>>> {
                static auto apply(std::pair<T*, std::vector<std::size_t>>& value) {
                    using alps::hdf5::get_pointer;
                    return value.first && (value.second.empty() || array_size(value.second)) ? get_pointer(*value.first) : nullptr;
                }
            };
            template<typename T> struct get_pointer<std::pair<T*, std::vector<std::size_t>> const> {
                static typename scalar_type<std::pair<T*, std::vector<std::size_t>>>::type const*
                apply(std::pair<T*, std::vector<std::size_t>> const& value) {
                    using alps::hdf5::get_pointer;
                    return value.first && (value.second.empty() || array_size(value.second)) ? get_pointer(*value.first) : nullptr;
                }
            };
        }

        template<typename T>
        void save(archive& ar, std::string const& path, std::pair<T*, std::vector<std::size_t>> const& value,
            std::vector<std::size_t> size = {}, std::vector<std::size_t> chunk = {}, std::vector<std::size_t> offset = {}) {
            if (is_continuous<T>::value) {
                auto extent = get_extent(value);
                size.insert(size.end(), extent.begin(), extent.end());
                chunk.insert(chunk.end(), extent.begin(), extent.end());
                offset.resize(offset.size() + extent.size(), 0);
                ar.write(path, get_pointer(value), size, chunk, offset);
            } else if (!detail::array_size(value.second)) {
                ar.write(path, static_cast<int const*>(nullptr), std::vector<std::size_t>());
            } else if (is_vectorizable(value)) {
                size.insert(size.end(), value.second.begin(), value.second.end());
                chunk.resize(chunk.size() + value.second.size(), 1);
                detail::visit_array(value.second, offset, [&](auto i, auto const& index) {
                    save(ar, path, value.first[i], size, chunk, index);
                });
            } else {
                if (path.find_last_of('@') != std::string::npos)
                    throw archive_error("attributes needs to be vectorizable: " + path + ALPS_STACKTRACE);
                if (ar.is_data(path)) ar.delete_data(path);
                detail::visit_array(value.second, {}, [&](auto i, auto const& index) {
                    save(ar, detail::array_path(path, index), value.first[i]);
                });
            }
        }

        template<typename T>
        void load(archive& ar, std::string const& path, std::pair<T*, std::vector<std::size_t>>& value,
            std::vector<std::size_t> chunk = {}, std::vector<std::size_t> offset = {}) {
            if (ar.is_group(path)) {
                detail::visit_array(value.second, {}, [&](auto i, auto const& index) {
                    load(ar, detail::array_path(path, index), value.first[i]);
                });
            } else {
                auto size = ar.extent(path);
                if (chunk.size() > size.size() || offset.size() > size.size())
                    throw archive_error("invalid data size" + ALPS_STACKTRACE);
                if (ar.is_null(path) && !detail::array_size(value.second)) return;
                set_extent(value, std::vector<std::size_t>(size.begin() + chunk.size(), size.end()));
                if (!value.second.empty() && !detail::array_size(value.second)) return;
                if (is_continuous<T>::value) {
                    chunk.insert(chunk.end(), size.begin() + chunk.size(), size.end());
                    offset.resize(size.size(), 0);
                    ar.read(path, get_pointer(value), chunk, offset);
                } else {
                    chunk.resize(chunk.size() + value.second.size(), 1);
                    detail::visit_array(value.second, offset, [&](auto i, auto const& index) {
                        load(ar, path, value.first[i], chunk, index);
                    });
                }
            }
        }
    }
}
#endif
