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

#ifndef ALPS_NGS_HDF5_HPP
#define ALPS_NGS_HDF5_HPP

#include <alps/ngs/config.hpp>
#include <alps/hdf5_export.h>
#include <alps/ngs/stacktrace.hpp>
#include <alps/hdf5/errors.hpp>
#include <alps/ngs/detail/remove_cvr.hpp>
#include <alps/ngs/detail/type_wrapper.hpp>

#include <boost/mpl/and.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/type_traits/is_array.hpp>
#include <boost/type_traits/remove_all_extents.hpp>


#include <complex>
#include <memory>
#include <vector>
#include <string>
#include <numeric>

#define ALPS_NGS_FOREACH_NATIVE_HDF5_TYPE(CALLBACK) \
    CALLBACK(char) \
    CALLBACK(signed char) \
    CALLBACK(unsigned char) \
    CALLBACK(short) \
    CALLBACK(unsigned short) \
    CALLBACK(int) \
    CALLBACK(unsigned) \
    CALLBACK(long) \
    CALLBACK(unsigned long) \
    CALLBACK(long long) \
    CALLBACK(unsigned long long) \
    CALLBACK(float) \
    CALLBACK(double) \
    CALLBACK(long double) \
    CALLBACK(bool) \
    CALLBACK(std::string) \
    CALLBACK(std::complex<float>) \
    CALLBACK(std::complex<double>) \
    CALLBACK(std::complex<long double>)

namespace alps {
    namespace hdf5 {

        namespace detail {
            struct archivecontext;
            class scoped_context;

            template<typename A, typename T> struct is_datatype_caller {
                static bool apply(A const & ar, std::string path) {
                    throw std::logic_error("only native datatypes can be probed: " + path + ALPS_STACKTRACE);
                    return false;
                }
            };

            #define ALPS_NGS_HDF5_IS_DATATYPE_CALLER(T) \
                template<typename A> struct is_datatype_caller<A, T > { \
                    static bool apply(A const & ar, std::string path, T unused = alps::detail::type_wrapper<T>::type()) { \
                        return ar.is_datatype_impl(path, unused); \
                    } \
                };
            ALPS_NGS_FOREACH_NATIVE_HDF5_TYPE(ALPS_NGS_HDF5_IS_DATATYPE_CALLER)
            #undef ALPS_NGS_HDF5_IS_DATATYPE_CALLER

            template<typename A> struct archive_proxy {

                explicit archive_proxy(std::string const & path, A & ar)
                    : path_(path), ar_(ar)
                {}

                template<typename T> archive_proxy & operator=(T const & value);
                template<typename T> archive_proxy & operator<<(T const & value);
                template <typename T> archive_proxy & operator>>(T & value);

                std::string path_;
                A ar_;
            };
        }

        class ALPS_HDF5_DECL archive {

            public:

                archive(boost::filesystem::path const & filename, std::string mode = "r");
                archive(archive const & arg);
                archive & operator=(archive const &) = delete;

                virtual ~archive();

                std::string const & get_filename() const;

                std::string encode_segment(std::string segment) const;
                std::string decode_segment(std::string segment) const;

                std::string get_context() const;
                void set_context(std::string const & context);
                std::string complete_path(std::string path) const;

                void close();
                bool is_open();

                bool is_data(std::string path) const;
                bool is_attribute(std::string path) const;
                bool is_group(std::string path) const;

                bool is_scalar(std::string path) const;
                bool is_null(std::string path) const;
                bool is_complex(std::string path) const;

                template<typename T> bool is_datatype(std::string path) const {
                    return detail::is_datatype_caller<archive, T>::apply(*this, path);
                }

                std::vector<std::string> list_children(std::string path) const;
                std::vector<std::string> list_attributes(std::string path) const;

                std::vector<std::size_t> extent(std::string path) const;
                std::size_t dimensions(std::string path) const;

                void create_group(std::string path) const;

                void delete_data(std::string path) const;
                void delete_group(std::string path) const;
                void delete_attribute(std::string path) const;


                detail::archive_proxy<archive> operator[](std::string const & path);

                template<typename T> void read(
                      std::string path
                    , T *
                    , std::vector<std::size_t>
                    , std::vector<std::size_t> = std::vector<std::size_t>()
                ) const {
                    throw std::logic_error("Invalid type on path: " + path + ALPS_STACKTRACE);
                }

                template<typename T> void write(
                      std::string path
                    , T const * value
                    , std::vector<std::size_t> size
                    , std::vector<std::size_t> chunk = std::vector<std::size_t>()
                    , std::vector<std::size_t> offset = std::vector<std::size_t>()
                ) const {
                    throw std::logic_error("Invalid type on path: " + path + ALPS_STACKTRACE);
                }

                #define ALPS_NGS_HDF5_DEFINE_API(T) \
                    void read(std::string path, T & value) const; \
                    void read( \
                          std::string path \
                        , T * value \
                        , std::vector<std::size_t> chunk \
                        , std::vector<std::size_t> offset = std::vector<std::size_t>() \
                    ) const; \
 \
                    void write(std::string path, T value) const; \
                    void write( \
                          std::string path \
                        , T const * value, std::vector<std::size_t> size \
                        , std::vector<std::size_t> chunk = std::vector<std::size_t>() \
                        , std::vector<std::size_t> offset = std::vector<std::size_t>() \
                    ) const;
                ALPS_NGS_FOREACH_NATIVE_HDF5_TYPE(ALPS_NGS_HDF5_DEFINE_API)
                #undef ALPS_NGS_HDF5_DEFINE_API

                #define ALPS_NGS_HDF5_IS_DATATYPE_IMPL_DECL(T) \
                    bool is_datatype_impl(std::string path, T) const;
                ALPS_NGS_FOREACH_NATIVE_HDF5_TYPE(ALPS_NGS_HDF5_IS_DATATYPE_IMPL_DECL)
                #undef ALPS_NGS_HDF5_IS_DATATYPE_IMPL_DECL

            private:

                friend class detail::scoped_context;
                std::string current_ = "/";
                std::shared_ptr<detail::archivecontext> context_;

        };

        template<typename T> struct is_continuous
            : public boost::false_type
        {};

        template<typename T> struct is_content_continuous
            : public is_continuous<T>
        {};


        template<typename T> struct scalar_type {
            typedef T type;
        };

        namespace detail {

            class scoped_context {
                public:
                    scoped_context(archive & ar, std::string const & path)
                        : ar_(ar), previous_(ar.get_context()) {
                        ar_.set_context(path);
                    }

                    ~scoped_context() noexcept {
                        // Restore without allocating while another exception
                        // may be unwinding from a user-defined save/load hook.
                        ar_.current_.swap(previous_);
                    }

                    scoped_context(scoped_context const &) = delete;
                    scoped_context & operator=(scoped_context const &) = delete;

                private:
                    archive & ar_;
                    std::string previous_;
            };

             template<typename T> struct get_extent {
                static std::vector<std::size_t> apply(T const & value) {
                    return std::vector<std::size_t>();
                }
            };

            template<typename T> struct set_extent {
                 static void apply(T &, std::vector<std::size_t> const &) {}
            };

            #define ALPS_NGS_HDF5_DEFINE_SET_EXTENT(T) \
                template<> struct set_extent<T> { \
                    static void apply(T &, std::vector<std::size_t> const & extent) { \
                        if (extent.size() > 0) \
                            throw wrong_type("The extents do not match" + ALPS_STACKTRACE); \
                    } \
                };
            ALPS_NGS_FOREACH_NATIVE_HDF5_TYPE(ALPS_NGS_HDF5_DEFINE_SET_EXTENT)
            #undef ALPS_NGS_HDF5_DEFINE_SET_EXTENT

            template<typename T> struct is_vectorizable {
                 static bool apply(T const & value){
                    return false;
                }
            };

            template<typename T> struct get_pointer {
                 static typename alps::hdf5::scalar_type<T>::type * apply(T &) {
                    return NULL;
                }
            };

            template<typename T> struct get_pointer<T const> {
                 static typename alps::hdf5::scalar_type<T>::type const * apply(T const &) {
                    return NULL;
                }
            };

        }

        template<typename T> typename scalar_type<T>::type * get_pointer(T & value) {
            return detail::get_pointer<T>::apply(value);
        }

        template<typename T> typename scalar_type<T>::type const * get_pointer(T const & value) {
            return detail::get_pointer<T const>::apply(value);
        }

        template<typename T> std::vector<std::size_t> get_extent(T const & value) {
            return detail::get_extent<T>::apply(value);
        }

        template<typename T> void set_extent(T & value, std::vector<std::size_t> const & size) {
            detail::set_extent<T>::apply(value, size);
        }

        template<typename T> bool is_vectorizable(T const & value) {
            return detail::is_vectorizable<T>::apply(value);
        }

        template<typename T> void save(
               archive & ar
             , std::string const & path
             , T const & value
             , std::vector<std::size_t> size = std::vector<std::size_t>()
             , std::vector<std::size_t> chunk = std::vector<std::size_t>()
             , std::vector<std::size_t> offset = std::vector<std::size_t>()
        ) {
            if (chunk.size())
                throw std::logic_error("user defined objects needs to be written continously" + ALPS_STACKTRACE);
            detail::scoped_context context(ar, path);
            value.save(ar);
        }

        template<typename T> void load(
               archive & ar
             , std::string const & path
             , T & value
             , std::vector<std::size_t> chunk = std::vector<std::size_t>()
             , std::vector<std::size_t> offset = std::vector<std::size_t>()
        ) {
            if (chunk.size())
                throw std::logic_error("user defined objects needs to be written continously" + ALPS_STACKTRACE);
            detail::scoped_context context(ar, path);
            value.load(ar);
        }

        #define ALPS_NGS_HDF5_DEFINE_FREE_FUNCTIONS(T) \
            template<> struct is_continuous< T > \
                : public boost::true_type \
            {}; \
            template<> struct is_continuous< T const > \
                : public boost::true_type \
            {}; \
 \
            namespace detail { \
                template<> struct ALPS_HDF5_DECL is_vectorizable< T > { \
                    static bool apply(T const & value); \
                }; \
                template<> struct ALPS_HDF5_DECL is_vectorizable< T const > { \
                    static bool apply(T & value); \
                }; \
 \
                template<> struct ALPS_HDF5_DECL get_pointer< T > { \
                    static alps::hdf5::scalar_type< T >::type * apply( T & value); \
                }; \
 \
                template<> struct ALPS_HDF5_DECL get_pointer< T const > { \
                    static alps::hdf5::scalar_type< T >::type const * apply( T const & value); \
                }; \
            } \
 \
            ALPS_HDF5_DECL void save( \
                  archive & ar \
                , std::string const & path \
                , T const & value \
                , std::vector<std::size_t> size = std::vector<std::size_t>() \
                , std::vector<std::size_t> chunk = std::vector<std::size_t>() \
                , std::vector<std::size_t> offset = std::vector<std::size_t>() \
            ); \
 \
            ALPS_HDF5_DECL void load( \
                  archive & ar \
                , std::string const & path \
                , T & value \
                , std::vector<std::size_t> chunk = std::vector<std::size_t>() \
                , std::vector<std::size_t> offset = std::vector<std::size_t>() \
            );
        ALPS_NGS_FOREACH_NATIVE_HDF5_TYPE(ALPS_NGS_HDF5_DEFINE_FREE_FUNCTIONS)
        #undef ALPS_NGS_HDF5_DEFINE_FREE_FUNCTIONS

        namespace detail {

            template<typename T> struct make_pvp_proxy {

                explicit make_pvp_proxy(std::string const & path, T value)
                    : path_(path), value_(value)
                {}

                make_pvp_proxy(make_pvp_proxy<T> const & arg)
                    : path_(arg.path_), value_(arg.value_)
                {}

                std::string path_;
                T value_;
            };

        }

        template <typename T> archive & operator<<(archive & ar, detail::make_pvp_proxy<T> const & proxy) {
            save(ar, proxy.path_, proxy.value_);
            return ar;
        }

        template <typename T> archive & operator>> (archive & ar, detail::make_pvp_proxy<T> proxy) {
            load(ar, proxy.path_, proxy.value_);
            return ar;
        }
    }

    template <typename T> typename boost::disable_if<typename boost::mpl::and_<
          typename boost::is_same<typename alps::detail::remove_cvr<typename boost::remove_all_extents<T>::type>::type, char>::type
        , typename boost::is_array<T>::type
    >::type, hdf5::detail::make_pvp_proxy<T &> >::type make_pvp(std::string const & path, T & value) {
        return hdf5::detail::make_pvp_proxy<T &>(path, value);
    }
    template <typename T> typename boost::disable_if<typename boost::mpl::and_<
          typename boost::is_same<typename alps::detail::remove_cvr<typename boost::remove_all_extents<T>::type>::type, char>::type
        , typename boost::is_array<T>::type
    >::type, hdf5::detail::make_pvp_proxy<T const &> >::type make_pvp(std::string const & path, T const & value) {
        return hdf5::detail::make_pvp_proxy<T const &>(path, value);
    }
    template <typename T> typename boost::enable_if<typename boost::mpl::and_<
          typename boost::is_same<typename alps::detail::remove_cvr<typename boost::remove_all_extents<T>::type>::type, char>::type
        , typename boost::is_array<T>::type
    >::type, hdf5::detail::make_pvp_proxy<std::string const> >::type make_pvp(std::string const & path, T const & value) {
        return hdf5::detail::make_pvp_proxy<std::string const>(path, value);
    }

    namespace hdf5 {
        namespace detail {

            template<typename A> template<typename T> archive_proxy<A> & archive_proxy<A>::operator=(T const & value) {
                ar_ << make_pvp(path_, value);
                return *this;
            }

            template<typename A> template<typename T> archive_proxy<A> & archive_proxy<A>::operator<<(T const & value) {
                return *this = value;
            }

            template<typename A> template <typename T> archive_proxy<A> & archive_proxy<A>::operator>> (T & value) {
                ar_ >> make_pvp(path_, value);
                return *this;
            }

        }
    }
}

#endif
