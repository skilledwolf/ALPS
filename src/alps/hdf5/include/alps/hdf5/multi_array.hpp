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

#ifndef ALPS_NGS_HDF5_BOOST_MULTI_ARRAY_HPP
#define ALPS_NGS_HDF5_BOOST_MULTI_ARRAY_HPP

#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/pair.hpp>

#include <alps/multi_array/multi_array.hpp>

#include <boost/multi_array.hpp>

namespace alps {
    namespace hdf5 {

        template<typename T, std::size_t N, typename A> struct scalar_type<boost::multi_array<T, N, A> > {
            typedef typename scalar_type<typename boost::remove_reference<typename boost::remove_cv<T>::type>::type>::type type;
        };
        template<typename T, std::size_t N, typename A> struct scalar_type<alps::multi_array<T, N, A> > 
            : public scalar_type<boost::multi_array<T, N, A> > 
        {};

        template<typename T, std::size_t N, typename A> struct is_content_continuous<alps::multi_array<T, N, A> > 
            : public is_continuous<T> 
        {};

        namespace detail {

            template<typename T, std::size_t N, typename A> struct get_extent<boost::multi_array<T, N, A> > {
                static std::vector<std::size_t> apply(boost::multi_array<T, N, A> const & value) {
                    using alps::hdf5::get_extent;
                    std::vector<std::size_t> result(value.shape(), value.shape() + boost::multi_array<T, N, A>::dimensionality);
                    if (value.num_elements()) {
                        std::vector<std::size_t> extent(get_extent(*value.data()));
                        for (std::size_t i = 1; i < value.num_elements(); ++i)
                            if (extent != get_extent(value.data()[i]))
                                throw archive_error("no rectengual matrix");
                        std::copy(extent.begin(), extent.end(), std::back_inserter(result));
                    }
                    return result;
                }
            };
            template<typename T, std::size_t N, typename A> struct get_extent<alps::multi_array<T, N, A> >
                : public get_extent<boost::multi_array<T, N, A> > 
            {};

            template<typename T, std::size_t N, typename A> struct set_extent<boost::multi_array<T, N, A> > {
                static void apply(boost::multi_array<T, N, A> & value, std::vector<std::size_t> const & size) {
                    using alps::hdf5::set_extent;
                    if ((is_continuous<T>::value && size.size() != N) || N > size.size())
                        throw archive_error("invalid data size");
                    if (!std::equal(value.shape(), value.shape() + boost::multi_array<T, N, A>::dimensionality, size.begin())) {
                        typename boost::multi_array<T, N, A>::extent_gen extents;
                        gen_extent(value, extents, size);
                    }
                    if (!is_continuous<T>::value && boost::multi_array<T, N, A>::dimensionality < size.size())
                        for (std::size_t i = 0; i < value.num_elements(); ++i)
                            set_extent(value.data()[i], std::vector<std::size_t>(size.begin() + boost::multi_array<T, N, A>::dimensionality, size.end()));
                }
                private:
                    template<std::size_t M> static void gen_extent(boost::multi_array<T, N, A> & value, boost::detail::multi_array::extent_gen<M> extents, std::vector<std::size_t> const & size) {
                        gen_extent(value, extents[size.front()], std::vector<std::size_t>(size.begin() + 1, size.end()));
                    }
                    static void gen_extent(boost::multi_array<T, N, A> & value, typename boost::detail::multi_array::extent_gen<N> extents, std::vector<std::size_t> const & size) {
                        value.resize(extents);
                    }
            };
            template<typename T, std::size_t N, typename A> struct set_extent<alps::multi_array<T, N, A> >
                : public set_extent<boost::multi_array<T, N, A> > 
            {};

            template<typename T, std::size_t N, typename A> struct is_vectorizable<boost::multi_array<T, N, A> > {
                static bool apply(boost::multi_array<T, N, A> const & value) {
                    using alps::hdf5::get_extent;
                    using alps::hdf5::is_vectorizable;
                    if (!value.num_elements())
                        return true;
                    std::vector<std::size_t> size(get_extent(*value.data()));
                    for (std::size_t i = 1; i < value.num_elements(); ++i)
                        if (!is_vectorizable(value.data()[i]) || size != get_extent(value.data()[i]))
                            return false;
                    return true;
                }
            };
            template<typename T, std::size_t N, typename A> struct is_vectorizable<alps::multi_array<T, N, A> >
                : public is_vectorizable<boost::multi_array<T, N, A> > 
            {};

            template<typename T, std::size_t N, typename A> struct get_pointer<boost::multi_array<T, N, A> > {
                static typename alps::hdf5::scalar_type<boost::multi_array<T, N, A> >::type * apply(boost::multi_array<T, N, A> & value) {
                    using alps::hdf5::get_pointer;
                    return value.num_elements() ? get_pointer(*value.data()) : nullptr;
                }
            };
            template<typename T, std::size_t N, typename A> struct get_pointer<alps::multi_array<T, N, A> >
                : public get_pointer<boost::multi_array<T, N, A> > 
            {};

            template<typename T, std::size_t N, typename A> struct get_pointer<boost::multi_array<T, N, A> const> {
                static typename alps::hdf5::scalar_type<boost::multi_array<T, N, A> >::type const * apply(boost::multi_array<T, N, A> const & value) {
                    using alps::hdf5::get_pointer;
                    return value.num_elements() ? get_pointer(*value.data()) : nullptr;
                }
            };
            template<typename T, std::size_t N, typename A> struct get_pointer<alps::multi_array<T, N, A> const>
                : public get_pointer<boost::multi_array<T, N, A> const> 
            {};

        }

        template<typename T, std::size_t N, typename A> void save(
              archive & ar
            , std::string const & path
            , boost::multi_array<T, N, A> const & value
            , std::vector<std::size_t> size = std::vector<std::size_t>()
            , std::vector<std::size_t> chunk = std::vector<std::size_t>()
            , std::vector<std::size_t> offset = std::vector<std::size_t>()
        ) {
            save(ar, path, std::make_pair(value.data(),
                 std::vector<std::size_t>(value.shape(), value.shape() + N)), size, chunk, offset);
        }
        template<typename T, std::size_t N, typename A> void save(
              archive & ar
            , std::string const & path
            , alps::multi_array<T, N, A> const & value
            , std::vector<std::size_t> size = std::vector<std::size_t>()
            , std::vector<std::size_t> chunk = std::vector<std::size_t>()
            , std::vector<std::size_t> offset = std::vector<std::size_t>()
        ) {
            save(ar, path, static_cast<boost::multi_array<T, N, A> const &>(value), size, chunk, offset);
        }

        template<typename T, std::size_t N, typename A> void load(
              archive & ar
            , std::string const & path
            , boost::multi_array<T, N, A> & value
            , std::vector<std::size_t> chunk = std::vector<std::size_t>()
            , std::vector<std::size_t> offset = std::vector<std::size_t>()
        ) {
            if (ar.is_group(path)) {
                if (!ar.list_children(path).empty())
                    throw invalid_path("invalid path");
                set_extent(value, std::vector<std::size_t>(N, 0));
            } else {
                const auto size = ar.extent(path);
                if (size.size() < chunk.size() + N)
                    throw archive_error("invalid data size");
                set_extent(value, std::vector<std::size_t>(size.begin() + chunk.size(), size.end()));
            }
            auto data = std::make_pair(value.data(),
                        std::vector<std::size_t>(value.shape(), value.shape() + N));
            load(ar, path, data, chunk, offset);
        }
        template<typename T, std::size_t N, typename A> void load(
              archive & ar
            , std::string const & path
            , alps::multi_array<T, N, A> & value
            , std::vector<std::size_t> chunk = std::vector<std::size_t>()
            , std::vector<std::size_t> offset = std::vector<std::size_t>()
        ) {
            load(ar, path, static_cast<boost::multi_array<T, N, A> &>(value), chunk, offset);                                                   
        }
   }
}

#endif
