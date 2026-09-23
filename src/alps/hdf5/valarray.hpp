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
#include <alps/hdf5/sequence.hpp>
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

        template<typename T> struct has_complex_elements<std::valarray<T> > 
            : public has_complex_elements<typename alps::detail::remove_cvr<T>::type>
        {};

        namespace detail {

            template<typename T> struct get_extent<std::valarray<T> > : sequence_get_extent<std::valarray<T>> {};

            template<typename T> struct set_extent<std::valarray<T> > : sequence_set_extent<std::valarray<T>> {};

            template<typename T> struct is_vectorizable<std::valarray<T> > : sequence_is_vectorizable<std::valarray<T>> {};

            template<typename T> struct get_pointer<std::valarray<T> > : sequence_get_pointer<std::valarray<T>> {};

            template<typename T> struct get_pointer<std::valarray<T> const> : sequence_get_pointer<std::valarray<T> const> {};
        }

        template<typename T> void save(
              archive & ar
            , std::string const & path
            , std::valarray<T> const & value
            , std::vector<std::size_t> size = std::vector<std::size_t>()
            , std::vector<std::size_t> chunk = std::vector<std::size_t>()
            , std::vector<std::size_t> offset = std::vector<std::size_t>()
        ) {
            detail::save_sequence(ar, path, value, std::move(size), std::move(chunk), std::move(offset));
        }

        template<typename T> void load(
              archive & ar
            , std::string const & path
            , std::valarray<T> & value
            , std::vector<std::size_t> chunk = std::vector<std::size_t>()
            , std::vector<std::size_t> offset = std::vector<std::size_t>()
        ) {
            detail::load_sequence(ar, path, value, std::move(chunk), std::move(offset));
        }
    }
}

#endif
