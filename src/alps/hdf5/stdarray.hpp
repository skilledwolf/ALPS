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

#ifndef ALPS_NGS_HDF5_STD_ARRAY_HPP
#define ALPS_NGS_HDF5_STD_ARRAY_HPP

#include <alps/hdf5/sequence.hpp>
#include <alps/ngs/cast.hpp>

#include <array>
#include <vector>
#include <iterator>
#include <algorithm>

namespace alps {
    namespace hdf5 {

        template<typename T, std::size_t N> struct scalar_type<std::array<T, N> > {
            typedef typename scalar_type<typename std::array<T, N>::value_type>::type type;
        };

        template<typename T, std::size_t N> struct is_continuous<std::array<T, N> >
            : public is_continuous<T>
        {};
        template<typename T, std::size_t N> struct is_continuous<std::array<T, N> const >
            : public is_continuous<T>
        {};

        template<typename T, std::size_t N> struct has_complex_elements<std::array<T, N> >
            : public has_complex_elements<typename alps::detail::remove_cvr<typename std::array<T, N>::value_type>::type>
        {};

        namespace detail {
            template<typename T, std::size_t N> struct get_extent<std::array<T, N>>
                : sequence_get_extent<std::array<T, N>, true> {};

            template<typename T, std::size_t N> struct set_extent<std::array<T, N>>
                : sequence_set_extent<std::array<T, N>, true> {};

            template<typename T, std::size_t N> struct is_vectorizable<std::array<T, N>>
                : sequence_is_vectorizable<std::array<T, N>, true> {};

            template<typename T, std::size_t N> struct get_pointer<std::array<T, N>>
                : sequence_get_pointer<std::array<T, N>> {};

            template<typename T, std::size_t N> struct get_pointer<std::array<T, N> const>
                : sequence_get_pointer<std::array<T, N> const> {};
        }

        template<typename T, std::size_t N> void save(
              archive & ar
            , std::string const & path
            , std::array<T, N> const & value
            , std::vector<std::size_t> size = std::vector<std::size_t>()
            , std::vector<std::size_t> chunk = std::vector<std::size_t>()
            , std::vector<std::size_t> offset = std::vector<std::size_t>()
        ) {
            detail::save_sequence(ar, path, value, std::move(size), std::move(chunk), std::move(offset));
        }

        template<typename T, std::size_t N> void load(
              archive & ar
            , std::string const & path
            , std::array<T, N> & value
            , std::vector<std::size_t> chunk = std::vector<std::size_t>()
            , std::vector<std::size_t> offset = std::vector<std::size_t>()
        ) {
            detail::load_sequence<std::array<T, N>, true>(ar, path, value, std::move(chunk), std::move(offset));
        }
    }
}

#endif
