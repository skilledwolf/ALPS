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

#ifndef ALPS_NGS_HDF5_BOOST_NUMERIC_UBLAS_MATRIX_HPP
#define ALPS_NGS_HDF5_BOOST_NUMERIC_UBLAS_MATRIX_HPP

#include <alps/hdf5/pair.hpp>

#include <boost/numeric/ublas/matrix.hpp>

#include <iterator>

namespace alps {
    namespace hdf5 {

        template <typename T, typename F, typename A> struct scalar_type<boost::numeric::ublas::matrix<T, F, A> > {
            typedef typename scalar_type<typename boost::remove_reference<typename boost::remove_cv<T>::type>::type>::type type;
        };

        template <typename T, typename F, typename A> struct has_complex_elements<boost::numeric::ublas::matrix<T, F, A> >
            : public has_complex_elements<typename alps::detail::remove_cvr<T>::type>
        {};

        namespace detail {
            // Preserve the historical physical-storage order for both uBLAS
            // layouts; the archive shape remains (rows, columns).
            template<class Matrix> auto matrix_view(Matrix& value) {
                return std::make_pair(value.size1() && value.size2() ? &value(0, 0) : nullptr,
                                      std::vector<std::size_t>{value.size1(), value.size2()});
            }

            template<typename T, typename F, typename A> struct get_extent<boost::numeric::ublas::matrix<T, F, A>> {
                static auto apply(boost::numeric::ublas::matrix<T, F, A> const& value) {
                    return alps::hdf5::get_extent(matrix_view(value));
                }
            };

            template<typename T, typename F, typename A> struct set_extent<boost::numeric::ublas::matrix<T, F, A>> {
                static void apply(boost::numeric::ublas::matrix<T, F, A>& value, std::vector<std::size_t> const& size) {
                    if (size.size() < 2) throw archive_error("invalid matrix dimensions" + ALPS_STACKTRACE);
                    value.resize(size[0], size[1], false);
                    auto view = matrix_view(value);
                    alps::hdf5::set_extent(view, size);
                }
            };

            template<typename T, typename F, typename A> struct is_vectorizable<boost::numeric::ublas::matrix<T, F, A>> {
                static bool apply(boost::numeric::ublas::matrix<T, F, A> const& value) {
                    return alps::hdf5::is_vectorizable(matrix_view(value));
                }
            };

            template<typename T, typename F, typename A> struct get_pointer<boost::numeric::ublas::matrix<T, F, A>> {
                static auto apply(boost::numeric::ublas::matrix<T, F, A>& value) {
                    auto view = matrix_view(value);
                    return alps::hdf5::get_pointer(view);
                }
            };

            template<typename T, typename F, typename A> struct get_pointer<boost::numeric::ublas::matrix<T, F, A> const> {
                static auto apply(boost::numeric::ublas::matrix<T, F, A> const& value) {
                    auto const view = matrix_view(value);
                    return alps::hdf5::get_pointer(view);
                }
            };
        }

        template <typename T, typename F, typename A> void save(
              archive & ar
            , std::string const & path
            , boost::numeric::ublas::matrix<T, F, A> const & value
            , std::vector<std::size_t> size = std::vector<std::size_t>()
            , std::vector<std::size_t> chunk = std::vector<std::size_t>()
            , std::vector<std::size_t> offset = std::vector<std::size_t>()
        ) {
            if constexpr (is_continuous<T>::value)
                save(ar, path, detail::matrix_view(value), size, chunk, offset);
            else
                throw wrong_type("invalid type" + ALPS_STACKTRACE);
        }

        template <typename T, typename F, typename A> void load(
              archive & ar
            , std::string const & path
            , boost::numeric::ublas::matrix<T, F, A> & value
            , std::vector<std::size_t> chunk = std::vector<std::size_t>()
            , std::vector<std::size_t> offset = std::vector<std::size_t>()
        ) {
            if constexpr (!is_continuous<T>::value) {
                throw invalid_path("invalid type" + ALPS_STACKTRACE);
            } else {
                if (ar.is_group(path)) throw invalid_path("invalid path" + ALPS_STACKTRACE);
                if (ar.is_null(path)) {
                    value.resize(0, 0, false);
                    return;
                }
                auto size = ar.extent(path);
                if (chunk.size() > size.size()) throw archive_error("invalid matrix dimensions" + ALPS_STACKTRACE);
                set_extent(value, std::vector<std::size_t>(size.begin() + chunk.size(), size.end()));
                auto view = detail::matrix_view(value);
                load(ar, path, view, chunk, offset);
            }
        }
    }
}

#endif
