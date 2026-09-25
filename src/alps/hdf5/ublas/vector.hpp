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

#ifndef ALPS_NGS_HDF5_BOOST_NUMERIC_UBLAS_VECTOR_HPP
#define ALPS_NGS_HDF5_BOOST_NUMERIC_UBLAS_VECTOR_HPP

#include <alps/hdf5/sequence.hpp>
#include <alps/ngs/cast.hpp>

#include <boost/numeric/ublas/vector.hpp>

#include <iterator>
#include <algorithm>

namespace alps {
    namespace hdf5 {

        template<typename T, typename A> struct scalar_type<boost::numeric::ublas::vector<T, A> > {
            typedef typename scalar_type<typename boost::numeric::ublas::vector<T, A>::value_type>::type type;
        };

        template<typename T, typename A> struct has_complex_elements<boost::numeric::ublas::vector<T, A> > 
            : public has_complex_elements<typename alps::detail::remove_cvr<typename boost::numeric::ublas::vector<T, A>::value_type>::type>
        {};

        namespace detail {

            template<typename T, typename A> struct get_extent<boost::numeric::ublas::vector<T, A> > : sequence_get_extent<boost::numeric::ublas::vector<T, A>> {};

            template<typename T, typename A> struct set_extent<boost::numeric::ublas::vector<T, A> > : sequence_set_extent<boost::numeric::ublas::vector<T, A>> {};

            template<typename T, typename A> struct is_vectorizable<boost::numeric::ublas::vector<T, A> > : sequence_is_vectorizable<boost::numeric::ublas::vector<T, A>> {};

            template<typename T, typename A> struct get_pointer<boost::numeric::ublas::vector<T, A> > : sequence_get_pointer<boost::numeric::ublas::vector<T, A>> {};

            template<typename T, typename A> struct get_pointer<boost::numeric::ublas::vector<T, A> const> : sequence_get_pointer<boost::numeric::ublas::vector<T, A> const> {};
        }

        template<typename T, typename A> void save(
              archive & ar
            , std::string const & path
            , boost::numeric::ublas::vector<T, A> const & value
            , std::vector<std::size_t> size = std::vector<std::size_t>()
            , std::vector<std::size_t> chunk = std::vector<std::size_t>()
            , std::vector<std::size_t> offset = std::vector<std::size_t>()
        ) {
            detail::save_sequence(ar, path, value, std::move(size), std::move(chunk), std::move(offset));
        }

        template<typename T, typename A> void load(
              archive & ar
            , std::string const & path
            , boost::numeric::ublas::vector<T, A> & value
            , std::vector<std::size_t> chunk = std::vector<std::size_t>()
            , std::vector<std::size_t> offset = std::vector<std::size_t>()
        ) {
            detail::load_sequence(ar, path, value, std::move(chunk), std::move(offset));
        }
    }
}

#endif
