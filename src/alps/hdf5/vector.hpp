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

#ifndef ALPS_NGS_HDF5_STD_VECTOR_HPP
#define ALPS_NGS_HDF5_STD_VECTOR_HPP

#include <alps/hdf5/sequence.hpp>
#include <alps/ngs/cast.hpp>

#include <boost/type_traits/is_scalar.hpp>

#include <vector>
#include <iterator>
#include <algorithm>

namespace alps {
    namespace hdf5 {

        template<typename T, typename A> struct scalar_type<std::vector<T, A> > {
            typedef typename scalar_type<typename std::vector<T, A>::value_type>::type type;
        };

        template<typename T, typename A> struct is_content_continuous<std::vector<T, A> >
            : public is_continuous<T> 
        {};
        template<typename A> struct is_content_continuous<std::vector<bool, A> >
            : public boost::false_type
        {};

        template<typename T, typename A> struct has_complex_elements<std::vector<T, A> > 
            : public has_complex_elements<typename alps::detail::remove_cvr<typename std::vector<T, A>::value_type>::type>
        {};

        namespace detail {

            template<typename T, typename A> struct get_extent<std::vector<T, A> > : sequence_get_extent<std::vector<T, A>> {};

            template<typename T, typename A> struct set_extent<std::vector<T, A> > : sequence_set_extent<std::vector<T, A>> {};

            template<typename A> struct set_extent<std::vector<bool, A> > {
                static void apply(std::vector<bool, A> & value, std::vector<std::size_t> const & extent) {
                    if (extent.size() != 1)
                        throw archive_error("dimensions do not match" + ALPS_STACKTRACE);
                    value.resize(extent[0]);
                }
            };

            template<typename T, typename A> struct is_vectorizable<std::vector<T, A> > : sequence_is_vectorizable<std::vector<T, A>> {};

            template<typename A> struct is_vectorizable<std::vector<bool, A> > {
                static bool apply(std::vector<bool, A> const & value) {
                    return true;
                }
            };

            template<typename T, typename A> struct get_pointer<std::vector<T, A> > : sequence_get_pointer<std::vector<T, A>> {};

            template<typename T, typename A> struct get_pointer<std::vector<T, A> const> : sequence_get_pointer<std::vector<T, A> const> {};

            template<typename A> struct get_pointer<std::vector<bool, A> > {
                static typename alps::hdf5::scalar_type<std::vector<bool, A> >::type * apply(std::vector<bool, A> & value) {
                    throw archive_error("std::vector<bool, A>[0] cannot be dereferenced" + ALPS_STACKTRACE);
                    return NULL;
                }
            };

            template<typename A> struct get_pointer<std::vector<bool, A> const> {
                static typename alps::hdf5::scalar_type<std::vector<bool, A> >::type const * apply(std::vector<bool, A> const & value) {
                    throw archive_error("std::vector<bool>[0] cannot be dereferenced" + ALPS_STACKTRACE);
                    return NULL;
                }
            };

        }


        template<typename T, typename A> void save(
              archive & ar
            , std::string const & path
            , std::vector<T, A> const & value
            , std::vector<std::size_t> size = std::vector<std::size_t>()
            , std::vector<std::size_t> chunk = std::vector<std::size_t>()
            , std::vector<std::size_t> offset = std::vector<std::size_t>()
        ) {
            detail::save_sequence(ar, path, value, std::move(size), std::move(chunk), std::move(offset));
        }

        template<typename A> void save(
              archive & ar
            , std::string const & path
            , std::vector<bool, A> const & value
            , std::vector<std::size_t> size = std::vector<std::size_t>()
            , std::vector<std::size_t> chunk = std::vector<std::size_t>()
            , std::vector<std::size_t> offset = std::vector<std::size_t>()
        ) {
            if (ar.is_group(path))
                ar.delete_group(path);
            if (value.size() == 0)
                ar.write(path, static_cast<bool const *>(NULL), std::vector<std::size_t>());
            else {
                size.push_back(value.size());
                chunk.push_back(1);
                offset.push_back(0);
                for(typename std::vector<bool, A>::const_iterator it = value.begin(); it != value.end(); ++it) {
                    offset.back() = it - value.begin();
                    bool const elem = *it;
                    ar.write(path, &elem, size, chunk, offset);
                }
            }
        }        

        template<typename T, typename A> void load(
              archive & ar
            , std::string const & path
            , std::vector<T, A> & value
            , std::vector<std::size_t> chunk = std::vector<std::size_t>()
            , std::vector<std::size_t> offset = std::vector<std::size_t>()
        ) {
            detail::load_sequence(ar, path, value, std::move(chunk), std::move(offset));
        }

        template<typename A> void load(
              archive & ar
            , std::string const & path
            , std::vector<bool, A> & value
            , std::vector<std::size_t> chunk = std::vector<std::size_t>()
            , std::vector<std::size_t> offset = std::vector<std::size_t>()
        ) {
            if (ar.is_group(path))
                throw archive_error("invalid dimensions" + ALPS_STACKTRACE);
            else {
                if (ar.is_complex(path))
                    throw archive_error("no complex value in archive" + ALPS_STACKTRACE);
                std::vector<std::size_t> size(ar.extent(path));
                if (size.size() == 0)
                    throw archive_error("invalid dimensions" + ALPS_STACKTRACE);
                else if (size[0] == 0)
                    value.resize(0);
                else {
                    value.resize(*(size.begin() + chunk.size()));
                    chunk.push_back(1);
                    offset.push_back(0);
                    for(typename std::vector<bool, A>::iterator it = value.begin(); it != value.end(); ++it) {
                        offset.back() = it - value.begin();
                        bool elem;
                        ar.read(path, &elem, chunk, offset);
                        *it = elem;
                    }
                }
            }
        }

    }
}

#endif
