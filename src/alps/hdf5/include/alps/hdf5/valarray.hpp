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

#ifndef ALPS_NGS_HDF5_STD_VALARRAY_HPP
#define ALPS_NGS_HDF5_STD_VALARRAY_HPP

#include <alps/hdf5.hpp>
#include <alps/ngs/cast.hpp>

#include <valarray>
#include <iterator>
#include <algorithm>

namespace alps {
    namespace hdf5 {

        template<typename T> struct scalar_type<std::valarray<T> > {
            typedef typename scalar_type<T>::type type;
        };

        template<typename T> struct is_content_continuous<std::valarray<T> >
            : public is_continuous<T> 
        {};

        namespace detail {

            template<typename T> struct get_extent<std::valarray<T> > {
                static std::vector<std::size_t> apply(std::valarray<T> const & value) {
                    using alps::hdf5::get_extent;
                    std::vector<std::size_t> result(1, value.size());
                    if (value.size()) {
                        std::vector<std::size_t> extent(get_extent(const_cast<std::valarray<T> &>(value)[0]));
                        if (!boost::is_scalar<T>::value)
                            for (std::size_t i = 1; i < value.size(); ++i)
                                if (extent != get_extent(const_cast<std::valarray<T> &>(value)[i]))
                                    throw archive_error("no rectengual matrix" + ALPS_STACKTRACE);
                        std::copy(extent.begin(), extent.end(), std::back_inserter(result));
                    } else if constexpr (is_continuous<T>::value) {
                        const auto element = get_extent(T{});
                        std::copy(element.begin(), element.end(), std::back_inserter(result));
                    }
                    return result;
                }
            };

            template<typename T> struct set_extent<std::valarray<T> > {
                static void apply(std::valarray<T> & value, std::vector<std::size_t> const & extent) {
                    using alps::hdf5::set_extent;
                    if (extent.empty())
                        throw archive_error("dimensions do not match" + ALPS_STACKTRACE);
                    value.resize(extent[0]);
                    if (extent.size() > 1)
                        for(std::size_t i = 0; i < value.size(); ++i)
                            set_extent(value[i], std::vector<std::size_t>(extent.begin() + 1, extent.end()));
                    else if (extent.size() == 0 && !boost::is_same<typename scalar_type<T>::type, T>::value)
                        throw archive_error("dimensions do not match" + ALPS_STACKTRACE);
                }
            };

            template<typename T> struct is_vectorizable<std::valarray<T> > {
                static bool apply(std::valarray<T> const & value) {
                    using alps::hdf5::get_extent;
                    using alps::hdf5::is_vectorizable;
                    if (value.size()) {
                        if (!is_vectorizable(const_cast<std::valarray<T> &>(value)[0]))
                            return false;
                        std::vector<std::size_t> first(get_extent(const_cast<std::valarray<T> &>(value)[0]));
                        if (!boost::is_scalar<T>::value) {
                            for(std::size_t i = 0; i < value.size(); ++i)
                                if (!is_vectorizable(const_cast<std::valarray<T> &>(value)[i])) {
                                    return false;
                                } else {
                                    std::vector<std::size_t> size(get_extent(const_cast<std::valarray<T> &>(value)[i]));
                                    if (
                                           first.size() != size.size() 
                                        || !std::equal(first.begin(), first.end(), size.begin())
                                    ) {
                                        return false;
                                    }
                                }
                        }
                    }
                    return true;
                }
            };

            template<typename T> struct get_pointer<std::valarray<T> > {
                static typename alps::hdf5::scalar_type<std::valarray<T> >::type * apply(std::valarray<T> & value) {
                    using alps::hdf5::get_pointer;
                    return value.size() ? get_pointer(value[0]) : nullptr;
                }
            };

            template<typename T> struct get_pointer<std::valarray<T> const> {
                static typename alps::hdf5::scalar_type<std::valarray<T> >::type const * apply(std::valarray<T> const & value) {
                    using alps::hdf5::get_pointer;
                    return value.size() ? get_pointer(const_cast<std::valarray<T> &>(value)[0]) : nullptr;
                }
            };
        }

        template<typename T> void save(
              archive & ar
            , std::string const & path
            , std::valarray<T> const & value
            , std::vector<std::size_t> size = std::vector<std::size_t>()
            , std::vector<std::size_t> chunk = std::vector<std::size_t>()
            , std::vector<std::size_t> offset = std::vector<std::size_t>()
        ) {
            if constexpr (is_continuous<T>::value) {
                std::vector<std::size_t> extent(get_extent(value));
                std::copy(extent.begin(), extent.end(), std::back_inserter(size));
                std::copy(extent.begin(), extent.end(), std::back_inserter(chunk));
                std::fill_n(std::back_inserter(offset), extent.size(), 0);
                ar.write(path, get_pointer(value), size, chunk, offset);
            } else if (value.size() == 0) {
                if (path.find_last_of('@') != std::string::npos)
                    throw archive_error("attributes need a native datatype" + ALPS_STACKTRACE);
                if (ar.is_group(path))
                    ar.delete_group(path);
                if (ar.is_data(path))
                    ar.delete_data(path);
                ar.create_group(path);
            }
            else if (is_vectorizable(value)) {
                size.push_back(value.size());
                chunk.push_back(1);
                offset.push_back(0);
                for(std::size_t i = 0; i < value.size(); ++i) {
                    offset.back() = i;
                    save(ar, path, const_cast<std::valarray<T> &>(value)[i], size, chunk, offset);
                }
            } else {
                if (path.find_last_of('@') != std::string::npos)
                    throw archive_error("attributes need a native datatype" + ALPS_STACKTRACE);
                if (ar.is_group(path))
                    ar.delete_group(path);
                if (ar.is_data(path))
                    ar.delete_data(path);
                ar.create_group(path);
                for(std::size_t i = 0; i < value.size(); ++i)
                    save(ar, ar.complete_path(path) + "/" + cast<std::string>(i), const_cast<std::valarray<T> &>(value)[i]);
            }
        }

        template<typename T> void load(
              archive & ar
            , std::string const & path
            , std::valarray<T> & value
            , std::vector<std::size_t> chunk = std::vector<std::size_t>()
            , std::vector<std::size_t> offset = std::vector<std::size_t>()
        ) {
            if (ar.is_group(path)) {
                std::vector<std::string> children = ar.list_children(path);
                for (const auto &child : children) {
                    const auto index = cast<std::size_t>(child);
                    if (index >= children.size() || cast<std::string>(index) != child)
                        throw invalid_path("invalid container index: " + child + ALPS_STACKTRACE);
                }
                value.resize(children.size());
                for (typename std::vector<std::string>::const_iterator it = children.begin(); it != children.end(); ++it)
                    load(ar, ar.complete_path(path) + "/" + *it, value[cast<std::size_t>(*it)]);
            } else {
                std::vector<std::size_t> size(ar.extent(path));
                if (size.size() <= chunk.size())
                    throw archive_error("dimensions do not match" + ALPS_STACKTRACE);
                if constexpr (is_continuous<T>::value) {
                    set_extent(value, std::vector<std::size_t>(size.begin() + chunk.size(), size.end()));
                    std::copy(size.begin() + chunk.size(), size.end(), std::back_inserter(chunk));
                    std::fill_n(std::back_inserter(offset), size.size() - offset.size(), 0);
                    ar.read(path, get_pointer(value), chunk, offset);
                } else {
                    set_extent(value, std::vector<std::size_t>(1, *(size.begin() + chunk.size())));
                    chunk.push_back(1);
                    offset.push_back(0);
                    for(std::size_t i = 0; i < value.size(); ++i) {
                        offset.back() = i;
                        load(ar, path, value[i], chunk, offset);
                    }
                }
            }
        }
    }
}

#endif
