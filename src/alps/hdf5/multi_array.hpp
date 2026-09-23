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

#include <alps/hdf5/pair.hpp>

#include <alps/multi_array.hpp>

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

        template<typename T, std::size_t N, typename A> struct has_complex_elements<boost::multi_array<T, N, A> > 
            : public has_complex_elements<typename alps::detail::remove_cvr<T>::type>
        {};
        template<typename T, std::size_t N, typename A> struct has_complex_elements<alps::multi_array<T, N, A> > 
            : public has_complex_elements<boost::multi_array<T, N, A> > 
        {};

        namespace detail {

            template<typename T, std::size_t N, typename A> struct get_extent<boost::multi_array<T, N, A> > {
                static std::vector<std::size_t> apply(boost::multi_array<T, N, A> const& value) {
                    auto view = std::make_pair(value.data(), std::vector<std::size_t>(value.shape(), value.shape() + N));
                    return alps::hdf5::get_extent(view);
                }
            };
            template<typename T, std::size_t N, typename A> struct get_extent<alps::multi_array<T, N, A> >
                : public get_extent<boost::multi_array<T, N, A> > 
            {};

            template<typename T, std::size_t N, typename A> struct set_extent<boost::multi_array<T, N, A> > {
                static void apply(boost::multi_array<T, N, A>& value, std::vector<std::size_t> const& size) {
                    if (size.size() < N) throw archive_error("invalid data size");
                    if (!std::equal(value.shape(), value.shape() + N, size.begin()))
                        value.resize(std::vector<std::size_t>(size.begin(), size.begin() + N));
                    auto view = std::make_pair(value.data(), std::vector<std::size_t>(value.shape(), value.shape() + N));
                    alps::hdf5::set_extent(view, size);
                }
            };
            template<typename T, std::size_t N, typename A> struct set_extent<alps::multi_array<T, N, A> >
                : public set_extent<boost::multi_array<T, N, A> > 
            {};

            template<typename T, std::size_t N, typename A> struct is_vectorizable<boost::multi_array<T, N, A> > {
                static bool apply(boost::multi_array<T, N, A> const& value) {
                    auto view = std::make_pair(value.data(), std::vector<std::size_t>(value.shape(), value.shape() + N));
                    return alps::hdf5::is_vectorizable(view);
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
            if (!is_continuous<T>::value && !is_vectorizable(value))
                throw wrong_type("invalid type");
            auto view = std::make_pair(value.data(), std::vector<std::size_t>(value.shape(), value.shape() + N));
            save(ar, path, view, size, chunk, offset);
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
            if (ar.is_group(path)) throw invalid_path("invalid path");
            if (ar.is_complex(path) != has_complex_elements<T>::value)
                throw archive_error("no complex value in archive" + ALPS_STACKTRACE);
            if (ar.is_null(path)) {
                value.resize(std::vector<std::size_t>(N, 0));
                return;
            }
            auto size = ar.extent(path);
            if (chunk.size() > size.size()) throw archive_error("invalid data size");
            set_extent(value, std::vector<std::size_t>(size.begin() + chunk.size(), size.end()));
            auto view = std::make_pair(value.data(), std::vector<std::size_t>(value.shape(), value.shape() + N));
            load(ar, path, view, chunk, offset);
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
