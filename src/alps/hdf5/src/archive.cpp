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

#include <alps/hdf5/archive.hpp>
#include <alps/ngs/cast.hpp>
#include <alps/ngs/config.hpp>
#include <alps/ngs/signal.hpp>
#include <alps/ngs/stacktrace.hpp>

#include <boost/scoped_array.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/noncopyable.hpp>

#include <hdf5.h>

#include <sstream>
#include <iostream>
#include <typeinfo>
#include <limits>
#include <type_traits>
#include <optional>

#ifdef ALPS_NGS_SINGLE_THREAD
    #define ALPS_HDF5_LOCK_MUTEX
#else
    #define ALPS_HDF5_LOCK_MUTEX boost::lock_guard<boost::recursive_mutex> guard(mutex_);
#endif

#ifdef H5_HAVE_THREADSAFE
    #define ALPS_HDF5_FAKE_THREADSAFETY
#else
    #define ALPS_HDF5_FAKE_THREADSAFETY ALPS_HDF5_LOCK_MUTEX
#endif

namespace alps {
    namespace hdf5 {

        namespace detail {

            template<typename T> struct native_ptr_converter {
                native_ptr_converter(std::size_t) {}
                inline T const * apply(T const * v) {
                    return v;
                }
            };

            template<> struct native_ptr_converter<std::string> {
                std::vector<char const *> data;
                native_ptr_converter(std::size_t size): data(size) {}
                inline char const * const * apply(std::string const * v) {
                    for (std::vector<char const *>::iterator it = data.begin(); it != data.end(); ++it)
                            *it = v[it - data.begin()].c_str();
                    return data.data();
                }
            };

            class ALPS_HDF5_DECL error {

                public:

                    std::string invoke(hid_t id) {
                        std::ostringstream buffer;
                        buffer << "HDF5 error: " << cast<std::string>(id) << std::endl;
                        H5Ewalk2(H5E_DEFAULT, H5E_WALK_DOWNWARD, callback, &buffer);
                        return buffer.str();
                    }

                private:

                    static herr_t callback(unsigned n, H5E_error2_t const * desc, void * buffer) {
                        *reinterpret_cast<std::ostringstream *>(buffer) 
                            << "    #" 
                            << cast<std::string>(n) 
                            << " " << desc->file_name 
                            << " line " 
                            << cast<std::string>(desc->line) 
                            << " in " 
                            << desc->func_name 
                            << "(): " 
                            << desc->desc 
                            << std::endl;
                        return 0;
                    }

            };

            template<herr_t(*F)(hid_t)> class resource : boost::noncopyable {
                public:
                    resource(hid_t id): _id(id) {
                        if (_id < 0)
                            throw archive_error(error().invoke(_id) + ALPS_STACKTRACE);
                    }

                    ~resource() {
                        if (_id >= 0 && F(_id) < 0) {
                            std::cerr << "Error in " 
                                      << __FILE__ 
                                      << " on " 
                                      << ALPS_NGS_STRINGIFY(__LINE__) 
                                      << " in " 
                                      << __FUNCTION__ // TODO: check for gcc and use __PRETTY_FUNCTION__ 
                                      << ":" 
                                      << std::endl 
                                      << error().invoke(_id) 
                                      << std::endl;
                            std::abort();
                        }
                    }

                    operator hid_t() const {
                        return _id;
                    }

                    hid_t release() noexcept {
                        hid_t id = _id;
                        _id = -1;
                        return id;
                    }

                private:
                    hid_t _id;
            };

            typedef resource<H5Gclose> group_type;
            typedef resource<H5Dclose> data_type;
            typedef resource<H5Aclose> attribute_type;
            typedef resource<H5Sclose> space_type;
            typedef resource<H5Tclose> type_type;
            typedef resource<H5Pclose> property_type;

            hid_t check_group(hid_t id) { group_type unused(id); return unused; }
            hid_t check_data(hid_t id) { data_type unused(id); return unused; }
            hid_t check_attribute(hid_t id) { attribute_type unused(id); return unused; }
            hid_t check_space(hid_t id) { space_type unused(id); return unused; }
            hid_t check_type(hid_t id) { type_type unused(id); return unused; }
            hid_t check_property(hid_t id) { property_type unused(id); return unused; }
            hid_t check_error(hid_t id) {
                if (id < 0)
                    throw archive_error(error().invoke(id) + ALPS_STACKTRACE);
                return id;
            }

            class vlen_cleanup : boost::noncopyable {
                public:
                    vlen_cleanup(hid_t type, hid_t space, void * buffer)
                        : type_(type), space_(space), buffer_(buffer) {}

                    ~vlen_cleanup() {
                        // Conversion may throw. Cleanup must not replace that
                        // exception with a reclamation error.
                        if (buffer_ != NULL)
                            H5Dvlen_reclaim(type_, space_, H5P_DEFAULT, buffer_);
                    }

                    void reclaim() {
                        void * buffer = buffer_;
                        buffer_ = NULL;
                        check_error(H5Dvlen_reclaim(type_, space_, H5P_DEFAULT, buffer));
                    }

                private:
                    hid_t type_, space_;
                    void * buffer_;
            };

            hid_t get_native_type(char) { return H5Tcopy(H5T_NATIVE_CHAR); }
            hid_t get_native_type(signed char) { return H5Tcopy(H5T_NATIVE_SCHAR); }
            hid_t get_native_type(unsigned char) { return H5Tcopy(H5T_NATIVE_UCHAR); }
            hid_t get_native_type(short) { return H5Tcopy(H5T_NATIVE_SHORT); }
            hid_t get_native_type(unsigned short) { return H5Tcopy(H5T_NATIVE_USHORT); }
            hid_t get_native_type(int) { return H5Tcopy(H5T_NATIVE_INT); }
            hid_t get_native_type(unsigned) { return H5Tcopy(H5T_NATIVE_UINT); }
            hid_t get_native_type(long) { return H5Tcopy(H5T_NATIVE_LONG); }
            hid_t get_native_type(unsigned long) { return H5Tcopy(H5T_NATIVE_ULONG); }
            hid_t get_native_type(long long) { return H5Tcopy(H5T_NATIVE_LLONG); }
            hid_t get_native_type(unsigned long long) { return H5Tcopy(H5T_NATIVE_ULLONG); }
            hid_t get_native_type(float) { return H5Tcopy(H5T_NATIVE_FLOAT); }
            hid_t get_native_type(double) { return H5Tcopy(H5T_NATIVE_DOUBLE); }
            hid_t get_native_type(long double) { return H5Tcopy(H5T_NATIVE_LDOUBLE); }
            hid_t get_native_type(bool) { return H5Tcopy(H5T_NATIVE_SCHAR); }
            hid_t get_native_type(std::string) {
                type_type type_id(H5Tcopy(H5T_C_S1));
                check_error(H5Tset_size(type_id, H5T_VARIABLE));
                return type_id.release();
            }

            // bool and signed char historically share H5T_NATIVE_SCHAR.
            // Mark new writes without changing the on-disk numeric type, so
            // old ALPS readers and native typed reads remain compatible.
            inline void mark_signed_byte(archive const & ar, std::string path,
                                         std::string const & kind) {
                std::size_t at = path.find_last_of('@');
                if (at != std::string::npos) {
                    std::string name = path.substr(at + 1);
                    if (name.compare(0, 2, "__") == 0)
                        return;
                    path = path.substr(0, at) + "@__alps_type__:" + name;
                } else
                    path += "/@__alps_type__";
                ar.write(path, kind);
            }
            template <typename T>
            void mark_signed_byte(archive const &, std::string const &, T *) {}
            inline void mark_signed_byte(archive const & ar, std::string const & path, bool *) {
                mark_signed_byte(ar, path, std::string("bool"));
            }
            inline void mark_signed_byte(archive const & ar, std::string const & path, signed char *) {
                mark_signed_byte(ar, path, std::string("int8"));
            }

            std::string attribute_parent(std::string const & path) {
                std::string parent = path.substr(0, path.find_last_of('@') - 1);
                return parent.empty() ? "/" : parent;
            }

            hid_t open_attribute(archive const & ar, hid_t file_id, std::string path) {
                if ((path = ar.complete_path(path)).find_last_of('@') == std::string::npos)
                    throw invalid_path("no attribute path: " + path + ALPS_STACKTRACE);
                return H5Aopen_by_name(file_id, detail::attribute_parent(path).c_str(), path.substr(path.find_last_of('@') + 1).c_str(), H5P_DEFAULT, H5P_DEFAULT);
            }

            herr_t list_children_visitor(hid_t, char const * n, const H5L_info_t *, void * d) {
                reinterpret_cast<std::vector<std::string> *>(d)->push_back(n);
                return 0;
            }

            herr_t list_attributes_visitor(hid_t, char const * n, const H5A_info_t *, void * d) {
                reinterpret_cast<std::vector<std::string> *>(d)->push_back(n);
                return 0;
            }

            struct ALPS_HDF5_DECL archivecontext : boost::noncopyable {

                archivecontext(std::string const & filename, bool write, bool replace, bool compress, bool large, bool memory)
                    : compress_(compress)
                    , write_(write || replace)
                    , replace_(!memory && replace)
                    , large_(large)
                    , memory_(memory)
                    , filename_(filename)
                {
                    construct();
                }

                ~archivecontext() {
                    destruct(true);
                }
                
                void grant(bool write, bool replace) {
                    if (!write_ && (write || replace)) {
                        destruct(false);
                        write_ = write || replace;
                        replace_ = !memory_ && replace;
                        construct();
                    }
                }
                
                bool compress_;
                bool write_;
                bool replace_;
                bool large_;
                bool memory_;
                std::string filename_;
                std::string suffix_;
                hid_t file_id_;
                
                private:

                    void construct() {
                        alps::ngs::signal::listen();
                        if (memory_ && large_)
                            throw archive_error("either memory or large file system can be used!" + ALPS_STACKTRACE);
                        else if (memory_) {
                            detail::property_type prop_id(H5Pcreate(H5P_FILE_ACCESS));
                            detail::check_error(H5Pset_fapl_core(prop_id, 1 << 20, true));
                            #ifndef ALPS_HDF5_CLOSE_GREEDY
                                detail::check_error(H5Pset_fclose_degree(prop_id, H5F_CLOSE_SEMI));
                            #endif
                            if (write_) {
                                if ((file_id_ = H5Fopen((filename_ + suffix_).c_str(), H5F_ACC_RDWR, prop_id)) < 0) {
                                    detail::property_type fcrt_id(H5Pcreate(H5P_FILE_CREATE));
                                    detail::check_error(H5Pset_link_creation_order(fcrt_id, (H5P_CRT_ORDER_TRACKED | H5P_CRT_ORDER_INDEXED)));
                                    detail::check_error(H5Pset_attr_creation_order(fcrt_id, (H5P_CRT_ORDER_TRACKED | H5P_CRT_ORDER_INDEXED)));
                                    detail::check_error(file_id_ = H5Fcreate((filename_ + suffix_).c_str(), H5F_ACC_TRUNC, fcrt_id, prop_id));
                                }
                            } else if ((file_id_ = H5Fopen((filename_ + suffix_).c_str(), H5F_ACC_RDONLY, prop_id)) < 0)
                                throw archive_not_found("file does not exists or is not a valid hdf5 archive: " + filename_ + suffix_ + ALPS_STACKTRACE);
                            else
                                detail::check_error(file_id_);
                        } else {
                            if (replace_ && large_)
                                throw archive_error("the combination 'wl' is not allowd!" + ALPS_STACKTRACE);
                            if (replace_)
                                for (std::size_t i = 0; boost::filesystem::exists(filename_ + (suffix_ = ".tmp." + cast<std::string>(i))); ++i);
                            if (write_ && replace_ && boost::filesystem::exists(filename_))
                                boost::filesystem::copy_file(filename_, filename_ + suffix_);
                            if (large_) {
                                {
                                    char filename0[4096], filename1[4096];
                                    snprintf(filename0, sizeof(filename0), filename_.c_str(), 0);
                                    snprintf(filename1, sizeof(filename1), filename_.c_str(), 1);
                                    if (!strcmp(filename0, filename1))
                                        throw archive_error("Large hdf5 archives need to have a '%d' part in the filename" + ALPS_STACKTRACE);
                                }
                                detail::property_type prop_id(H5Pcreate(H5P_FILE_ACCESS));
                                detail::check_error(H5Pset_fapl_family(prop_id, 1 << 30, H5P_DEFAULT));
                                #ifndef ALPS_HDF5_CLOSE_GREEDY
                                    detail::check_error(H5Pset_fclose_degree(prop_id, H5F_CLOSE_SEMI));
                                #endif
                                if (write_) {
                                    if ((file_id_ = H5Fopen((filename_ + suffix_).c_str(), H5F_ACC_RDWR, prop_id)) < 0) {
                                        detail::property_type fcrt_id(H5Pcreate(H5P_FILE_CREATE));
                                        detail::check_error(H5Pset_link_creation_order(fcrt_id, (H5P_CRT_ORDER_TRACKED | H5P_CRT_ORDER_INDEXED)));
                                        detail::check_error(H5Pset_attr_creation_order(fcrt_id, (H5P_CRT_ORDER_TRACKED | H5P_CRT_ORDER_INDEXED)));
                                        detail::check_error(file_id_ = H5Fcreate((filename_ + suffix_).c_str(), H5F_ACC_TRUNC, fcrt_id, prop_id));
                                    }
                                } else if ((file_id_ = H5Fopen((filename_ + suffix_).c_str(), H5F_ACC_RDONLY, prop_id)) < 0)
                                    throw archive_not_found("file does not exists or is not a valid hdf5 archive: " + filename_ + suffix_ + ALPS_STACKTRACE);
                                else
                                    detail::check_error(file_id_);
                            } else {
                                if (!write_) {
                                    if (!boost::filesystem::exists(filename_ + suffix_))
                                        throw archive_not_found("file does not exist: " + filename_ + suffix_ + ALPS_STACKTRACE);
                                    if (detail::check_error(H5Fis_hdf5((filename_ + suffix_).c_str())) == 0)
                                        throw archive_error("no valid hdf5 file: " + filename_ + suffix_ + ALPS_STACKTRACE);
                                }
                                #ifndef ALPS_HDF5_CLOSE_GREEDY
                                    detail::property_type ALPS_HDF5_FILE_ACCESS(H5Pcreate(H5P_FILE_ACCESS));
                                    detail::check_error(H5Pset_fclose_degree(ALPS_HDF5_FILE_ACCESS, H5F_CLOSE_SEMI));
                                #else
                                    #define ALPS_HDF5_FILE_ACCESS H5P_DEFAULT
                                #endif
                                if (write_) {
                                    if ((file_id_ = H5Fopen((filename_ + suffix_).c_str(), H5F_ACC_RDWR, ALPS_HDF5_FILE_ACCESS)) < 0) {
                                        detail::property_type fcrt_id(H5Pcreate(H5P_FILE_CREATE));
                                        detail::check_error(H5Pset_link_creation_order(fcrt_id, (H5P_CRT_ORDER_TRACKED | H5P_CRT_ORDER_INDEXED)));
                                        detail::check_error(H5Pset_attr_creation_order(fcrt_id, (H5P_CRT_ORDER_TRACKED | H5P_CRT_ORDER_INDEXED)));
                                        detail::check_error(file_id_ = H5Fcreate((filename_ + suffix_).c_str(), H5F_ACC_TRUNC, fcrt_id, ALPS_HDF5_FILE_ACCESS));
                                    }
                                } else
                                    detail::check_error(file_id_ = H5Fopen((filename_ + suffix_).c_str(), H5F_ACC_RDONLY, ALPS_HDF5_FILE_ACCESS));
                                #ifdef ALPS_HDF5_CLOSE_GREEDY
                                    #undef(ALPS_HDF5_FILE_ACCESS)
                                #endif
                            }
                        }
                    }

                    void destruct(bool abort) {
                        try {
                            H5Fflush(file_id_, H5F_SCOPE_GLOBAL);
                            #ifndef ALPS_HDF5_CLOSE_GREEDY
                                if (
                                       H5Fget_obj_count(file_id_, H5F_OBJ_DATATYPE) > 0
                                    || H5Fget_obj_count(file_id_, H5F_OBJ_ALL) - H5Fget_obj_count(file_id_, H5F_OBJ_FILE) > 0
                                ) {
                                    std::cerr << "Not all resources closed in file '" << filename_ << suffix_ << "'" << std::endl;
                                    std::abort();
                                }
                            #endif
                            if (H5Fclose(file_id_) < 0)
                                std::cerr << "Error in " 
                                          << __FILE__ 
                                          << " on " 
                                          << ALPS_NGS_STRINGIFY(__LINE__) 
                                          << " in " 
                                          << __FUNCTION__ // TODO: check for gcc and use __PRETTY_FUNCTION__  
                                          << ":" 
                                          << std::endl
                                          << error().invoke(file_id_)
                                          << std::endl;
                            if (replace_) {
                                if (boost::filesystem::exists(filename_))
                                    boost::filesystem::remove(filename_);
                                boost::filesystem::rename(filename_ + suffix_, filename_);
                            }
                        } catch (std::exception & ex) {
                            if (abort) {
                                std::cerr << "Error destroying HDF5 context of file '" << filename_ << suffix_ << "'\n" << ex.what() << std::endl;
                                std::abort();
                            } else
                                throw ex;
                        }
                    }
            };

        }

        archive::archive(std::string const & filename, int props) { // TODO: remove that!
            construct(filename, props);
        }

        archive::archive(std::string const & filename, char prop) { // TODO: remove that!
            construct(filename,
                  ('w' == prop ? WRITE | REPLACE : 0)
                | ('a' == prop ? WRITE : 0)
                | ('c' == prop ? COMPRESS : 0)
                | ('l' == prop ? LARGE : 0)
                | ('m' == prop ? MEMORY : 0)
            );
        }

        archive::archive(std::string const & filename, char signed prop) { // TODO: remove that!
            construct(filename,
                  ('w' == prop ? WRITE | REPLACE : 0)
                | ('a' == prop ? WRITE : 0)
                | ('c' == prop ? COMPRESS : 0)
                | ('l' == prop ? LARGE : 0)
                | ('m' == prop ? MEMORY : 0)
            );
        }

        archive::archive(boost::filesystem::path const & filename, std::string mode) {
            construct(filename.string(),
                  (mode.find_last_of('w') == std::string::npos ? 0 : WRITE | REPLACE)
                | (mode.find_last_of('a') == std::string::npos ? 0 : WRITE)
                | (mode.find_last_of('c') == std::string::npos ? 0 : COMPRESS)
                | (mode.find_last_of('l') == std::string::npos ? 0 : LARGE)
                | (mode.find_last_of('m') == std::string::npos ? 0 : MEMORY)
            );
            
        }

        archive::archive(archive const & arg)
            : current_(arg.current_)
            , context_(arg.context_)
        {
            if (context_ != NULL) {
                ALPS_HDF5_LOCK_MUTEX
                ++ref_cnt_[file_key(context_->filename_, context_->large_, context_->memory_)].second;
            }
        }

        archive::~archive() {
            if (context_ != NULL)
                try {
                    close();
                } catch (std::exception & ex) {
                    std::cerr << "Error destructing archive of file '" << ex.what() << std::endl;
                    std::abort();
                }
        }

        void archive::abort() {
            // Do not use a lock here, else deadlocking is really likly
            for (std::map<std::string, std::pair<detail::archivecontext *, std::size_t> >::iterator it = ref_cnt_.begin(); it != ref_cnt_.end(); ++it) {
                bool replace = it->second.first->replace_;
                std::string filename = it->second.first->filename_;
                it->second.first->replace_ = false;
                delete it->second.first;
                if (replace && boost::filesystem::exists(filename))
                    boost::filesystem::remove(filename);
            }
            ref_cnt_.clear();
        }

        void archive::close() {
            if (context_ == NULL)
                throw archive_closed("the archive is closed" + ALPS_STACKTRACE);
            ALPS_HDF5_LOCK_MUTEX
            H5Fflush(context_->file_id_, H5F_SCOPE_GLOBAL);
            if (!--ref_cnt_[file_key(context_->filename_, context_->large_, context_->memory_)].second) {
                ref_cnt_.erase(file_key(context_->filename_, context_->large_, context_->memory_));
                delete context_;
            }
            context_ = NULL;
        }
        
        bool archive::is_open() {
            return context_ != NULL;
        }

        std::string const & archive::get_filename() const {
            if (context_ == NULL)
                throw archive_closed("the archive is closed" + ALPS_STACKTRACE);
            return context_->filename_;
        }

        std::string archive::encode_segment(std::string segment) const {
            char chars[] = {'&', '/'};
            for (std::size_t i = 0; i < sizeof(chars); ++i)
                for (std::size_t pos = segment.find_first_of(chars[i]); pos < std::string::npos; pos = segment.find_first_of(chars[i], pos + 1))
                    segment = segment.substr(0, pos) + "&#" + cast<std::string>(static_cast<int>(chars[i])) + ";" + segment.substr(pos + 1);
            return segment;
        }

        std::string archive::decode_segment(std::string segment) const {
            for (std::size_t pos = segment.find_first_of('&'); pos < std::string::npos; pos = segment.find_first_of('&', pos + 1))
                segment = segment.substr(0, pos) 
                        + static_cast<char>(cast<int>(segment.substr(pos + 2, segment.find_first_of(';', pos) - pos - 2))) 
                        + segment.substr(segment.find_first_of(';', pos) + 1);
            return segment;
        }

        std::string archive::get_context() const {
            return current_;
        }
    
        void archive::set_context(std::string const & context) {
            ALPS_HDF5_LOCK_MUTEX
            current_ = complete_path(context);
        }
    
        std::string archive::complete_path(std::string path) const {
            if (path.size() > 1 && *path.rbegin() == '/')
                path = path.substr(0, path.size() - 1);
            if (path.size() && path[0] == '/')
                return path;
            else if (path.size() < 2 || path.substr(0, 2) != "..")
                return current_ + (current_.size() == 1 || !path.size() ? "" : "/") + path;
            else {
                std::string ctx = current_;
                while (ctx.size() && path.size() && path.substr(0, 2) == "..") {
                    ctx = ctx.substr(0, ctx.find_last_of('/'));
                    path = path.size() == 2 ? "" : path.substr(3);
                }
                if (ctx.empty()) ctx = "/";
                return ctx + (ctx.size() == 1 || !path.size() ? "" : "/") + path;
            }
        }
    
        bool archive::is_data(std::string path) const {
            if (context_ == NULL)
                throw archive_closed("the archive is closed" + ALPS_STACKTRACE);
            if ((path = complete_path(path)).find_last_of('@') != std::string::npos)
                throw invalid_path("no data path: " + path + ALPS_STACKTRACE);
            ALPS_HDF5_FAKE_THREADSAFETY
            hid_t id = H5Dopen2(context_->file_id_, path.c_str(), H5P_DEFAULT);
            return id < 0 ? false : detail::check_data(id) != 0;
        }
    
        bool archive::is_attribute(std::string path) const {
            if (context_ == NULL)
                throw archive_closed("the archive is closed" + ALPS_STACKTRACE);
            if ((path = complete_path(path)).find_last_of('@') == std::string::npos)
                return false;
            ALPS_HDF5_FAKE_THREADSAFETY
            return detail::check_error(H5Aexists_by_name(context_->file_id_, detail::attribute_parent(path).c_str(), path.substr(path.find_last_of('@') + 1).c_str(), H5P_DEFAULT));
        }
    
        bool archive::is_group(std::string path) const {
            if (context_ == NULL)
                throw archive_closed("the archive is closed" + ALPS_STACKTRACE);
            if ((path = complete_path(path)).find_last_of('@') != std::string::npos)
                return false;
            ALPS_HDF5_FAKE_THREADSAFETY
            hid_t id = H5Gopen2(context_->file_id_, path.c_str(), H5P_DEFAULT);
            return id < 0 ? false : detail::check_group(id) != 0;
        }
    
        bool archive::is_scalar(std::string path) const {
            hid_t space_id;
            if (context_ == NULL)
                throw archive_closed("the archive is closed" + ALPS_STACKTRACE);
            ALPS_HDF5_FAKE_THREADSAFETY
            if ((path = complete_path(path)).find_last_of('@') != std::string::npos && is_attribute(path)) {
                detail::attribute_type attr_id(detail::open_attribute(*this, context_->file_id_, path));
                space_id = H5Aget_space(attr_id);
            } else if (path.find_last_of('@') == std::string::npos && is_data(path)) {
                detail::data_type data_id(H5Dopen2(context_->file_id_, path.c_str(), H5P_DEFAULT));
                space_id = H5Dget_space(data_id);
            } else
                #ifdef ALPS_HDF5_READ_GREEDY
                    return false;
                #else
                    throw path_not_found("error reading path " + path + ALPS_STACKTRACE);
                #endif
            H5S_class_t type = H5Sget_simple_extent_type(space_id);
            detail::check_space(space_id);
            if (type == H5S_NO_CLASS)
                throw archive_error("error reading class " + path + ALPS_STACKTRACE);
            return type == H5S_SCALAR;
        }

        bool archive::is_null(std::string path) const {
            hid_t space_id;
            if (context_ == NULL)
                throw archive_closed("the archive is closed" + ALPS_STACKTRACE);
            ALPS_HDF5_FAKE_THREADSAFETY
            if ((path = complete_path(path)).find_last_of('@') != std::string::npos) {
                detail::attribute_type attr_id(detail::open_attribute(*this, context_->file_id_, path));
                space_id = H5Aget_space(attr_id);
            } else {
                detail::data_type data_id(H5Dopen2(context_->file_id_, path.c_str(), H5P_DEFAULT));
                space_id = H5Dget_space(data_id);
            }
            H5S_class_t type = H5Sget_simple_extent_type(space_id);
            detail::check_space(space_id);
            if (type == H5S_NO_CLASS)
                throw archive_error("error reading class " + path + ALPS_STACKTRACE);
            return type == H5S_NULL;
        }
    
        bool archive::is_complex(std::string path) const {
            if (context_ == NULL)
                throw archive_closed("the archive is closed" + ALPS_STACKTRACE);
            ALPS_HDF5_FAKE_THREADSAFETY
            if ((path = complete_path(path)).find_last_of('@') != std::string::npos)
                return is_attribute(path.substr(0, path.find_last_of('@')) + "@__complex__:" + path.substr(path.find_last_of('@') + 1))
                    && is_scalar(path.substr(0, path.find_last_of('@')) + "@__complex__:" + path.substr(path.find_last_of('@') + 1));
            else if (is_group(path)) {
                std::vector<std::string> children = list_children(path);
                for (std::size_t i = 0; i < children.size(); ++i)
                    if (is_complex(path + "/" + children[i]))
                        return true;
                return false;
            } else
                return is_attribute(path + "/@__complex__") && is_scalar(path + "/@__complex__");
        }
    
        std::vector<std::string> archive::list_children(std::string path) const {
            if (context_ == NULL)
                throw archive_closed("the archive is closed" + ALPS_STACKTRACE);
            if ((path = complete_path(path)).find_last_of('@') != std::string::npos)
                throw invalid_path("no group path: " + path + ALPS_STACKTRACE);
            std::vector<std::string> list;
            ALPS_HDF5_FAKE_THREADSAFETY
            if (!is_group(path))
                throw path_not_found("The group '" + path + "' does not exist." + ALPS_STACKTRACE);
            detail::group_type group_id(H5Gopen2(context_->file_id_, path.c_str(), H5P_DEFAULT));
            detail::check_error(H5Literate(group_id, H5_INDEX_NAME, H5_ITER_NATIVE, NULL, detail::list_children_visitor, &list));
            return list;
        }
    
        std::vector<std::string> archive::list_attributes(std::string path) const {
            if (context_ == NULL)
                throw archive_closed("the archive is closed" + ALPS_STACKTRACE);
            if ((path = complete_path(path)).find_last_of('@') != std::string::npos)
                throw invalid_path("no group or data path: " + path + ALPS_STACKTRACE);
            std::vector<std::string> list;
            ALPS_HDF5_FAKE_THREADSAFETY
            if (is_group(path)) {
                detail::group_type id(H5Gopen2(context_->file_id_, path.c_str(), H5P_DEFAULT));
                detail::check_error(H5Aiterate2(id, H5_INDEX_CRT_ORDER, H5_ITER_NATIVE, NULL, detail::list_attributes_visitor, &list));
            } else if (is_data(path)) {
                detail::data_type id(H5Dopen2(context_->file_id_, path.c_str(), H5P_DEFAULT));
                detail::check_error(H5Aiterate2(id, H5_INDEX_CRT_ORDER, H5_ITER_NATIVE, NULL, detail::list_attributes_visitor, &list));
            } else
                throw path_not_found("The path '" + path + "' does not exist." + ALPS_STACKTRACE);
            return list;
        }
    
        std::vector<std::size_t> archive::extent(std::string path) const {
            if (context_ == NULL)
                throw archive_closed("the archive is closed" + ALPS_STACKTRACE);
            if (is_null(path = complete_path(path)))
                return std::vector<std::size_t>(1, 0);
            else if (is_scalar(path))
                return std::vector<std::size_t>(1, 1);
            std::vector<hsize_t> buffer(dimensions(path), 0);
            hid_t space_id;
            ALPS_HDF5_FAKE_THREADSAFETY
            if (path.find_last_of('@') != std::string::npos) {
                detail::attribute_type attr_id(detail::open_attribute(*this, context_->file_id_, path));
                space_id = H5Aget_space(attr_id);
            } else {
                detail::data_type data_id(H5Dopen2(context_->file_id_, path.c_str(), H5P_DEFAULT));
                space_id = H5Dget_space(data_id);
            }
            detail::check_error(H5Sget_simple_extent_dims(space_id, &buffer.front(), NULL));
            detail::check_space(space_id);
            std::vector<std::size_t> extent(buffer.begin(), buffer.end());
            return extent;
        }
    
        std::size_t archive::dimensions(std::string path) const {
            if (context_ == NULL)
                throw archive_closed("the archive is closed" + ALPS_STACKTRACE);
            ALPS_HDF5_FAKE_THREADSAFETY
            if ((path = complete_path(path)).find_last_of('@') != std::string::npos) {
                detail::attribute_type attr_id(detail::open_attribute(*this, context_->file_id_, path));
                return detail::check_error(H5Sget_simple_extent_dims(detail::space_type(H5Aget_space(attr_id)), NULL, NULL));
            } else {
                detail::data_type data_id(H5Dopen2(context_->file_id_, path.c_str(), H5P_DEFAULT));
                return detail::check_error(H5Sget_simple_extent_dims(detail::space_type(H5Dget_space(data_id)), NULL, NULL));
            }
        }
    
        void archive::create_group(std::string path) const {
            if (context_ == NULL)
                throw archive_closed("the archive is closed" + ALPS_STACKTRACE);
            if ((path = complete_path(path)).find_last_of('@') != std::string::npos)
                throw invalid_path("no group path: " + path + ALPS_STACKTRACE);
            ALPS_HDF5_FAKE_THREADSAFETY
            if (is_data(path))
                delete_data(path);
            if (!is_group(path)) {
                std::size_t pos;
                hid_t group_id = -1;
                for (pos = path.find_last_of('/'); group_id < 0 && pos > 0 && pos < std::string::npos; pos = path.find_last_of('/', pos - 1))
                    group_id = H5Gopen2(context_->file_id_, path.substr(0, pos).c_str(), H5P_DEFAULT);
                if (group_id < 0) {
                    if ((pos = path.find_first_of('/', 1)) != std::string::npos) {
                        detail::property_type prop_id(H5Pcreate(H5P_GROUP_CREATE));
                        detail::check_error(H5Pset_link_creation_order(prop_id, (H5P_CRT_ORDER_TRACKED | H5P_CRT_ORDER_INDEXED)));
                        detail::check_error(H5Pset_attr_creation_order(prop_id, (H5P_CRT_ORDER_TRACKED | H5P_CRT_ORDER_INDEXED)));
                        detail::check_group(H5Gcreate2(context_->file_id_, path.substr(0, pos).c_str(), H5P_DEFAULT, prop_id, H5P_DEFAULT));
                    }
                } else {
                    pos = path.find_first_of('/', pos + 1);
                    detail::check_group(group_id);
                }
                while (pos != std::string::npos && (pos = path.find_first_of('/', pos + 1)) != std::string::npos && pos > 0) {
                    detail::property_type prop_id(H5Pcreate(H5P_GROUP_CREATE));
                    detail::check_error(H5Pset_link_creation_order(prop_id, (H5P_CRT_ORDER_TRACKED | H5P_CRT_ORDER_INDEXED)));
                    detail::check_error(H5Pset_attr_creation_order(prop_id, (H5P_CRT_ORDER_TRACKED | H5P_CRT_ORDER_INDEXED)));
                    detail::check_group(H5Gcreate2(context_->file_id_, path.substr(0, pos).c_str(), H5P_DEFAULT, prop_id, H5P_DEFAULT));
                }                
                detail::property_type prop_id(H5Pcreate(H5P_GROUP_CREATE));
                detail::check_error(H5Pset_link_creation_order(prop_id, (H5P_CRT_ORDER_TRACKED | H5P_CRT_ORDER_INDEXED)));
                detail::check_error(H5Pset_attr_creation_order(prop_id, (H5P_CRT_ORDER_TRACKED | H5P_CRT_ORDER_INDEXED)));
                detail::check_group(H5Gcreate2(context_->file_id_, path.c_str(), H5P_DEFAULT, prop_id, H5P_DEFAULT));
            }
        }
    
        void archive::delete_data(std::string path) const {
            if (context_ == NULL)
                throw archive_closed("the archive is closed" + ALPS_STACKTRACE);
            if ((path = complete_path(path)).find_last_of('@') != std::string::npos)
                throw invalid_path("no data path: " + path + ALPS_STACKTRACE);
            ALPS_HDF5_FAKE_THREADSAFETY
            if (is_data(path))
                detail::check_error(H5Ldelete(context_->file_id_, path.c_str(), H5P_DEFAULT));
            else if (is_group(path))
                throw invalid_path("the path contains a group: " + path + ALPS_STACKTRACE);
        }
    
        void archive::delete_group(std::string path) const  {
            if (context_ == NULL)
                throw archive_closed("the archive is closed" + ALPS_STACKTRACE);
            if ((path = complete_path(path)).find_last_of('@') != std::string::npos)
                throw invalid_path("no group path: " + path + ALPS_STACKTRACE);
            ALPS_HDF5_FAKE_THREADSAFETY
            if (is_group(path))
                detail::check_error(H5Ldelete(context_->file_id_, path.c_str(), H5P_DEFAULT));
            else if (is_data(path))
                throw invalid_path("the path contains a dataset: " + path + ALPS_STACKTRACE);
        }
    
        void archive::delete_attribute(std::string path) const {
            if (context_ == NULL)
                throw archive_closed("the archive is closed" + ALPS_STACKTRACE);
            if ((path = complete_path(path)).find_last_of('@') == std::string::npos)
                throw invalid_path("no attribute path: " + path + ALPS_STACKTRACE);
            ALPS_HDF5_FAKE_THREADSAFETY
            if (is_attribute(path)) {
                auto at = path.find_last_of('@');
                detail::check_error(H5Adelete_by_name(context_->file_id_,
                    detail::attribute_parent(path).c_str(), path.substr(at + 1).c_str(), H5P_DEFAULT));
            }
        }
    
        void archive::set_complex(std::string path) {
            if (context_ == NULL)
                throw archive_closed("the archive is closed" + ALPS_STACKTRACE);
            ALPS_HDF5_FAKE_THREADSAFETY
            // Resolve against the current context first, as every sibling here
            // does. Without it a relative path -- in particular the empty path
            // that paramvalue's saver uses, `ar[""] << value` under a context
            // set to the parameter name -- produced the marker attribute path
            // "/@__complex__", i.e. an attribute on the root group rather than
            // on the dataset just written. Saving a complex parameter then
            // failed with "HDF5 error: -1".
            path = complete_path(path);
            if (path.find_last_of('@') != std::string::npos)
                write(path.substr(0, path.find_last_of('@')) + "@__complex__:" + path.substr(path.find_last_of('@') + 1), true);
            else {
                if (is_group(path)) {
                    std::vector<std::string> children = list_children(path);
                    for (std::vector<std::string>::const_iterator it = children.begin(); it != children.end(); ++it)
                        set_complex(path + "/" + *it);
                } else
                    write(path + "/@__complex__", true);
            }
        }
    
        detail::archive_proxy<archive> archive::operator[](std::string const & path) {
            return detail::archive_proxy<archive>(path, *this);
        }
    
        namespace detail {

            // A read owns one dataset or attribute, its types and its file space.
            // Only opening, closing and transferring bytes depend on the object kind.
            herr_t close_storage_object(hid_t id) {
                return H5Iget_type(id) == H5I_ATTR ? H5Aclose(id) : H5Dclose(id);
            }

            using storage_type = resource<close_storage_object>;

            struct read_source : boost::noncopyable {
                read_source(archive const & ar, hid_t file, std::string const & path)
                    : attribute(path.find_last_of('@') != std::string::npos)
                    , object(attribute ? open_attribute(ar, file, path)
                                       : H5Dopen2(file, path.c_str(), H5P_DEFAULT))
                    , type(attribute ? H5Aget_type(object) : H5Dget_type(object))
                    , native_type(H5Tget_native_type(type, H5T_DIR_ASCEND))
                    , space(attribute ? H5Aget_space(object) : H5Dget_space(object))
                {}

                void read(hid_t memory_space, void * buffer) const {
                    check_error(attribute ? H5Aread(object, native_type, buffer)
                                          : H5Dread(object, native_type, memory_space,
                                                    space, H5P_DEFAULT, buffer));
                }

                bool const attribute;
                storage_type object;
                type_type type, native_type;
                space_type space;
            };

            template <bool Scalar, typename T>
            void read_buffer(read_source const & source, hid_t memory_space,
                             T * value, std::size_t count, std::string const & path) {
                if (H5Tget_class(source.native_type) == H5T_STRING) {
                    if (!check_error(H5Tis_variable_str(source.type))) {
                        if constexpr (Scalar) {
                            std::string raw(H5Tget_size(source.type) + 1, '\0');
                            source.read(memory_space, &raw[0]);
                            *value = cast<T>(raw);
                        } else
                            throw std::logic_error("multidimensional dataset of fixed string datas is not implemented (" + path + ")" + ALPS_STACKTRACE);
                    } else {
                        boost::scoped_array<char *> raw(new char *[count]);
                        source.read(memory_space, raw.get());
                        vlen_cleanup cleanup(source.native_type, memory_space, raw.get());
                        // Preserve the existing conversion contract: scalar strings
                        // are parsed as std::string; array elements use char* casts.
                        if constexpr (Scalar)
                            *value = cast<T>(std::string(raw[0]));
                        else
                            cast(raw.get(), raw.get() + count, value);
                        cleanup.reclaim();
                    }
                    return;
                }
                // Dispatch once on the stored numeric type for every read shape.
                // Conversion remains in ALPS rather than HDF5 to preserve cast semantics.
                #define ALPS_HDF5_READ_NUMERIC(U)                          \
                    if (check_error(H5Tequal(source.native_type, type_type(get_native_type(static_cast<U>(0)))))) { \
                        boost::scoped_array<U> raw(new U[count]);                  \
                        source.read(memory_space, raw.get());                     \
                        cast(raw.get(), raw.get() + count, value);                 \
                        return;                                                   \
                    }
                ALPS_HDF5_READ_NUMERIC(char)
                ALPS_HDF5_READ_NUMERIC(signed char)
                ALPS_HDF5_READ_NUMERIC(unsigned char)
                ALPS_HDF5_READ_NUMERIC(short)
                ALPS_HDF5_READ_NUMERIC(unsigned short)
                ALPS_HDF5_READ_NUMERIC(int)
                ALPS_HDF5_READ_NUMERIC(unsigned)
                ALPS_HDF5_READ_NUMERIC(long)
                ALPS_HDF5_READ_NUMERIC(unsigned long)
                ALPS_HDF5_READ_NUMERIC(long long)
                ALPS_HDF5_READ_NUMERIC(unsigned long long)
                ALPS_HDF5_READ_NUMERIC(float)
                ALPS_HDF5_READ_NUMERIC(double)
                ALPS_HDF5_READ_NUMERIC(long double)
                ALPS_HDF5_READ_NUMERIC(bool)
                #undef ALPS_HDF5_READ_NUMERIC
                throw wrong_type("invalid type" + ALPS_STACKTRACE);
            }

            template <typename T>
            void read_scalar(archive const & ar, hid_t file, std::string const & path, T & value) {
                bool const attribute = path.find_last_of('@') != std::string::npos;
                if (!(attribute ? ar.is_attribute(path) : ar.is_data(path)))
                    throw path_not_found("the path does not exist: " + path + ALPS_STACKTRACE);
                if (!ar.is_scalar(path))
                    throw wrong_type("scalar - vector conflict in path: " + path + ALPS_STACKTRACE);
                read_source source(ar, file, path);
                read_buffer<true>(source, source.space, &value, 1, path);
            }

            std::size_t checked_product(std::vector<std::size_t> const & dimensions,
                                        std::size_t initial = 1) {
                if (std::find(dimensions.begin(), dimensions.end(), 0) != dimensions.end())
                    return 0;
                for (auto dimension : dimensions) {
                    if (dimension > (std::numeric_limits<std::size_t>::max)() / initial)
                        throw archive_error("HDF5 extent product overflows" + ALPS_STACKTRACE);
                    initial *= dimension;
                }
                return initial;
            }

            void validate_selection(std::vector<std::size_t> const & size,
                                    std::vector<std::size_t> const & chunk,
                                    std::vector<std::size_t> const & offset,
                                    std::string const & path) {
                if (size.size() != chunk.size() || size.size() != offset.size())
                    throw archive_error("wrong chunk or offset passed for path: " + path + ALPS_STACKTRACE);
                for (std::size_t i = 0; i < size.size(); ++i)
                    if (offset[i] > size[i] || chunk[i] > size[i] - offset[i])
                        throw archive_error("selection exceeds data size for path: " + path + ALPS_STACKTRACE);
            }

            template <typename T>
            void read_array(archive const & ar, hid_t file, std::string const & path, T * value,
                            std::vector<std::size_t> const & chunk, std::vector<std::size_t> offset) {
                std::vector<std::size_t> const size = ar.extent(path);
                if (offset.empty())
                    offset.resize(ar.dimensions(path), 0);
                validate_selection(size, chunk, offset, path);
                if (ar.is_null(path))
                    return;
                for (auto extent : chunk)
                    if (extent == 0)
                        throw archive_error("size is zero in one dimension in path: " + path + ALPS_STACKTRACE);
                bool const attribute = path.find_last_of('@') != std::string::npos;
                if (!(attribute ? ar.is_attribute(path) : ar.is_data(path)))
                    throw path_not_found("the path does not exist: " + path + ALPS_STACKTRACE);
                if (ar.is_scalar(path)) {
                    if (attribute)
                        throw wrong_type("scalar - vector conflict in path: " + path + ALPS_STACKTRACE);
                    throw archive_error("scalar - vector conflict in path: " + path + ALPS_STACKTRACE);
                }
                read_source source(ar, file, path);
                if (attribute && chunk != size) {
                    if (H5Tget_class(source.native_type) == H5T_STRING)
                        throw std::logic_error("non continous multidimensional dataset as attributes are not implemented (" + path + ")" + ALPS_STACKTRACE);
                    throw std::logic_error("Not Implemented, path: " + path + ALPS_STACKTRACE);
                }
                std::vector<hsize_t> const dimensions(chunk.begin(), chunk.end());
                space_type memory_space(H5Screate_simple(static_cast<int>(dimensions.size()), dimensions.data(), NULL));
                if (!attribute) {
                    std::vector<hsize_t> const start(offset.begin(), offset.end());
                    check_error(H5Sselect_hyperslab(source.space, H5S_SELECT_SET, start.data(), NULL, dimensions.data(), NULL));
                }
                std::size_t const count = checked_product(chunk);
                read_buffer<false>(source, memory_space, value, count, path);
            }
        }

        // Keep the exported overloads; all read logic lives in ordinary functions.
        #define ALPS_HDF5_DEFINE_READ(T)                                          \
            void archive::read(std::string path, T & value) const {               \
                ALPS_HDF5_FAKE_THREADSAFETY                                       \
                if (context_ == NULL)                                            \
                    throw archive_closed("the archive is closed" + ALPS_STACKTRACE); \
                detail::read_scalar(*this, context_->file_id_, complete_path(path), value); \
            }                                                                    \
            void archive::read(std::string path, T * value, std::vector<std::size_t> chunk, std::vector<std::size_t> offset) const { \
                ALPS_HDF5_FAKE_THREADSAFETY                                       \
                if (context_ == NULL)                                            \
                    throw archive_closed("the archive is closed" + ALPS_STACKTRACE); \
                detail::read_array(*this, context_->file_id_, complete_path(path), value, chunk, offset); \
            }
        ALPS_NGS_FOREACH_NATIVE_HDF5_TYPE(ALPS_HDF5_DEFINE_READ)
        #undef ALPS_HDF5_DEFINE_READ

        namespace detail {
            void dataset_layout(hid_t properties, std::vector<std::size_t> const & size,
                                std::size_t element_size, bool compress) {
                auto const bytes = checked_product(size, element_size);
                check_error(H5Pset_fill_time(properties, H5D_FILL_TIME_NEVER));
                if (bytes < ALPS_HDF5_SZIP_BLOCK_SIZE * element_size)
                    check_error(H5Pset_layout(properties, H5D_COMPACT));
                else {
                    // HDF5 chunks must be smaller than 4 GiB. Halving the largest
                    // dimension also works when leading dimensions are singletons.
                    auto chunk = size;
                    while (checked_product(chunk, element_size) >= (1ULL << 32)) {
                        auto largest = std::max_element(chunk.begin(), chunk.end());
                        *largest = *largest / 2 + *largest % 2;
                    }
                    std::vector<hsize_t> const dimensions(chunk.begin(), chunk.end());
                    check_error(H5Pset_chunk(properties, static_cast<int>(dimensions.size()), dimensions.data()));
                }
                if (compress && bytes > ALPS_HDF5_SZIP_BLOCK_SIZE * element_size)
                    check_error(H5Pset_szip(properties, H5_SZIP_NN_OPTION_MASK, ALPS_HDF5_SZIP_BLOCK_SIZE));
            }

            template <typename T>
            hid_t open_write_target(archive const & ar, hid_t location, std::string const & path,
                                    hid_t type, hid_t space, std::vector<std::size_t> const & size,
                                    bool compress) {
                auto const at = path.find_last_of('@');
                bool const attribute = at != std::string::npos;
                std::string const parent = attribute ? attribute_parent(path) : path.substr(0, path.find_last_of('/'));
                std::string const name = attribute ? path.substr(at + 1) : path;
                if (!attribute) {
                    if (ar.is_group(path))
                        ar.delete_group(path);
                    if (!parent.empty())
                        ar.create_group(parent);
                }
                hid_t const existing = attribute
                    ? H5Aopen(location, name.c_str(), H5P_DEFAULT)
                    : H5Dopen2(location, path.c_str(), H5P_DEFAULT);
                if (existing >= 0) {
                    storage_type current(existing);
                    space_type current_space(attribute ? H5Aget_space(current) : H5Dget_space(current));
                    auto const kind = H5Sget_simple_extent_type(space);
                    bool same_shape = H5Sget_simple_extent_type(current_space) == kind;
                    if (same_shape && kind == H5S_SIMPLE) {
                        auto rank = check_error(H5Sget_simple_extent_ndims(current_space));
                        std::vector<hsize_t> dimensions(static_cast<std::size_t>(rank));
                        check_error(H5Sget_simple_extent_dims(current_space, dimensions.data(), NULL));
                        same_shape = dimensions.size() == size.size()
                            && std::equal(dimensions.begin(), dimensions.end(), size.begin());
                    }
                    type_type stored_type(attribute ? H5Aget_type(current) : H5Dget_type(current));
                    type_type native_type(H5Tget_native_type(stored_type, H5T_DIR_ASCEND));
                    bool const same_type = std::is_same<T, std::string>::value
                        ? H5Tget_class(native_type) == H5T_STRING
                        : check_error(H5Tequal(native_type, type)) > 0;
                    if (same_shape && same_type)
                        return current.release();
                }
                if (existing >= 0) {
                    if (attribute)
                        check_error(H5Adelete(location, name.c_str()));
                    else
                        check_error(H5Ldelete(location, path.c_str(), H5P_DEFAULT));
                }
                if (attribute)
                    return H5Acreate2(location, name.c_str(), type, space, H5P_DEFAULT, H5P_DEFAULT);
                property_type properties(H5Pcreate(H5P_DATASET_CREATE));
                check_error(H5Pset_attr_creation_order(properties, H5P_CRT_ORDER_TRACKED | H5P_CRT_ORDER_INDEXED));
                if (H5Sget_simple_extent_type(space) == H5S_SIMPLE && !std::is_same<T, std::string>::value)
                    dataset_layout(properties, size, sizeof(T), compress);
                return H5Dcreate2(location, path.c_str(), type, space, H5P_DEFAULT, properties, H5P_DEFAULT);
            }

            template <typename T>
            void write_buffer(archive const & ar, archivecontext const & context,
                              std::string const & path, T const * value, bool scalar,
                              std::vector<std::size_t> const & size = {},
                              std::vector<std::size_t> chunk = {}, std::vector<std::size_t> offset = {}) {
                bool const attribute = path.find_last_of('@') != std::string::npos;
                bool const null = !scalar && std::all_of(size.begin(), size.end(), [](auto n) { return n == 0; });
                std::size_t count = 1;
                if (!scalar) {
                    if (chunk.empty()) chunk = size;
                    if (offset.empty()) offset.resize(size.size(), 0);
                    validate_selection(size, chunk, offset, path);
                    checked_product(size, sizeof(T));
                    count = null ? 0 : checked_product(chunk);
                    if (attribute && chunk != size)
                        throw std::logic_error("partial attribute writes are not supported: " + path + ALPS_STACKTRACE);
                }
                if (count && !value)
                    throw archive_error("null write buffer for path: " + path + ALPS_STACKTRACE);
                std::vector<hsize_t> const dimensions(size.begin(), size.end());
                space_type space(scalar ? H5Screate(H5S_SCALAR) : null ? H5Screate(H5S_NULL)
                    : H5Screate_simple(static_cast<int>(dimensions.size()), dimensions.data(), NULL));
                type_type type(get_native_type(T()));
                native_ptr_converter<T> converter(count);
                // Keep the parent object alive through writes to existing
                // attributes as well as creation and nested type-marker writes.
                std::optional<resource<H5Oclose>> parent;
                if (attribute) {
                    auto const parent_path = attribute_parent(path);
                    if (!ar.is_group(parent_path) && !ar.is_data(parent_path))
                        throw path_not_found("unknown path: " + parent_path + ALPS_STACKTRACE);
                    parent.emplace(H5Oopen(context.file_id_, parent_path.c_str(), H5P_DEFAULT));
                }
                storage_type target(open_write_target<T>(ar, attribute ? static_cast<hid_t>(*parent) : context.file_id_,
                                                        path, type, space, size, context.compress_));
                if (count) {
                    if (attribute)
                        check_error(H5Awrite(target, type, converter.apply(value)));
                    else if (scalar)
                        check_error(H5Dwrite(target, type, space, space, H5P_DEFAULT, converter.apply(value)));
                    else {
                        std::vector<hsize_t> const start(offset.begin(), offset.end()), extent(chunk.begin(), chunk.end());
                        check_error(H5Sselect_hyperslab(space, H5S_SELECT_SET, start.data(), NULL, extent.data(), NULL));
                        space_type memory(H5Screate_simple(static_cast<int>(extent.size()), extent.data(), NULL));
                        check_error(H5Dwrite(target, type, memory, space, H5P_DEFAULT, converter.apply(value)));
                    }
                }
                mark_signed_byte(ar, path, static_cast<T *>(nullptr));
            }
        }

        #define ALPS_HDF5_DEFINE_WRITE(T)                                         \
            void archive::write(std::string path, T value) const {                \
                ALPS_HDF5_FAKE_THREADSAFETY                                        \
                if (context_ == NULL)                                             \
                    throw archive_closed("the archive is closed" + ALPS_STACKTRACE); \
                if (!context_->write_)                                           \
                    throw archive_error("the archive is not writeable" + ALPS_STACKTRACE); \
                detail::write_buffer(*this, *context_, complete_path(path), &value, true); \
            }                                                                     \
            void archive::write(std::string path, T const * value, std::vector<std::size_t> size, std::vector<std::size_t> chunk, std::vector<std::size_t> offset) const { \
                ALPS_HDF5_FAKE_THREADSAFETY                                        \
                if (context_ == NULL)                                             \
                    throw archive_closed("the archive is closed" + ALPS_STACKTRACE); \
                if (!context_->write_)                                           \
                    throw archive_error("the archive is not writeable" + ALPS_STACKTRACE); \
                detail::write_buffer(*this, *context_, complete_path(path), value, false, size, chunk, offset); \
            }
        ALPS_NGS_FOREACH_NATIVE_HDF5_TYPE(ALPS_HDF5_DEFINE_WRITE)
        #undef ALPS_HDF5_DEFINE_WRITE

        #define ALPS_NGS_HDF5_IMPLEMENT_FREE_FUNCTIONS(T)                                                                                                               \
            namespace detail {                                                                                                                                          \
                alps::hdf5::scalar_type< T >::type * get_pointer< T >::apply( T & value) {                                                                              \
                    return &value;                                                                                                                                      \
                }                                                                                                                                                       \
                                                                                                                                                                        \
                alps::hdf5::scalar_type< T >::type const * get_pointer< T const >::apply( T const & value) {                                                            \
                    return &value;                                                                                                                                      \
                }                                                                                                                                                       \
                                                                                                                                                                        \
                bool is_vectorizable< T >::apply(T const & value) {                                                                                                     \
                    return true;                                                                                                                                        \
                }                                                                                                                                                       \
                bool is_vectorizable< T const >::apply(T & value) {                                                                                                     \
                    return true;                                                                                                                                        \
                }                                                                                                                                                       \
            }                                                                                                                                                           \
                                                                                                                                                                        \
            void save(                                                                                                                                                  \
                    archive & ar                                                                                                                                        \
                , std::string const & path                                                                                                                              \
                , T const & value                                                                                                                                       \
                , std::vector<std::size_t> size                                                                                                                         \
                , std::vector<std::size_t> chunk                                                                                                                        \
                , std::vector<std::size_t> offset                                                                                                                       \
            ){                                                                                                                                                          \
                if (!size.size())                                                                                                                                       \
                    ar.write(path, value);                                                                                                                              \
                else                                                                                                                                                    \
                    ar.write(path, get_pointer(value), size, chunk, offset);                                                                                            \
            }                                                                                                                                                           \
                                                                                                                                                                        \
            void load(                                                                                                                                                  \
                  archive & ar                                                                                                                                          \
                , std::string const & path                                                                                                                              \
                , T & value                                                                                                                                             \
                , std::vector<std::size_t> chunk                                                                                                                        \
                , std::vector<std::size_t> offset                                                                                                                       \
            ) {                                                                                                                                                         \
                if (!chunk.size())                                                                                                                                      \
                    ar.read(path, value);                                                                                                                               \
                else                                                                                                                                                    \
                    ar.read(path, get_pointer(value), chunk, offset);                                                                                                   \
            }
        ALPS_NGS_FOREACH_NATIVE_HDF5_TYPE(ALPS_NGS_HDF5_IMPLEMENT_FREE_FUNCTIONS)
        #undef ALPS_NGS_HDF5_IMPLEMENT_FREE_FUNCTIONS

        namespace detail {
            template<typename T> struct is_datatype_impl_compare {
                static bool apply(type_type const & native_id) {
                    return check_error(
                        H5Tequal(type_type(H5Tcopy(native_id)), type_type(get_native_type(typename alps::detail::type_wrapper<T>::type())))
                    ) > 0;
                }
            };
            template<> struct is_datatype_impl_compare<std::string> {
                static bool apply(type_type const & native_id) {
                    return H5Tget_class(native_id) == H5T_STRING;
                }
            };
        }

        #define ALPS_NGS_HDF5_IS_DATATYPE_IMPL_IMPL(T)                                                                                                                  \
            bool archive::is_datatype_impl(std::string path, T) const {                                                                                                 \
                ALPS_HDF5_FAKE_THREADSAFETY                                                                                                                             \
                hid_t type_id;                                                                                                                                          \
                path = complete_path(path);                                                                                                                             \
                if (context_ == NULL)                                                                                                                                   \
                    throw archive_closed("the archive is closed" + ALPS_STACKTRACE);                                                                                    \
                if (path.find_last_of('@') != std::string::npos && is_attribute(path)) {                                                                                \
                    detail::attribute_type attr_id(detail::open_attribute(*this, context_->file_id_, path));                                                            \
                    type_id = H5Aget_type(attr_id);                                                                                                                     \
                } else if (path.find_last_of('@') == std::string::npos && is_data(path)) {                                                                              \
                    detail::data_type data_id(H5Dopen2(context_->file_id_, path.c_str(), H5P_DEFAULT));                                                                 \
                    type_id = H5Dget_type(data_id);                                                                                                                     \
                } else                                                                                                                                                  \
                    throw path_not_found("no valid path: " + path + ALPS_STACKTRACE);                                                                                   \
                detail::type_type native_id(H5Tget_native_type(type_id, H5T_DIR_ASCEND));                                                                               \
                detail::check_type(type_id);                                                                                                                            \
                {                                                                                                                                                       \
                    ALPS_HDF5_LOCK_MUTEX                                                                                                                                \
                    return detail::is_datatype_impl_compare< T >::apply(native_id);                                                                                     \
                }                                                                                                                                                       \
            }
        ALPS_NGS_FOREACH_NATIVE_HDF5_TYPE(ALPS_NGS_HDF5_IS_DATATYPE_IMPL_IMPL)
        #undef ALPS_NGS_HDF5_IS_DATATYPE_IMPL_IMPL

        void archive::construct(std::string const & filename, std::size_t props) {
            // Contexts must be absolute so restoring the root after nested
            // save/load calls cannot resolve an empty string relative to the
            // last dataset visited.
            current_ = "/";
            ALPS_HDF5_LOCK_MUTEX
            detail::check_error(H5Eset_auto2(H5E_DEFAULT, NULL, NULL));
            if (props & COMPRESS) {
                unsigned int flag;
                detail::check_error(H5Zget_filter_info(H5Z_FILTER_SZIP, &flag));
                props &= (flag & H5Z_FILTER_CONFIG_ENCODE_ENABLED ? ~0x00 : ~COMPRESS);
            }
            if (ref_cnt_.find(file_key(filename, props & LARGE, props & MEMORY)) == ref_cnt_.end())
                ref_cnt_.insert(std::make_pair(
                      file_key(filename, props & LARGE, props & MEMORY)
                    , std::make_pair(context_ = new detail::archivecontext(filename, props & WRITE, props & REPLACE, props & COMPRESS, props & LARGE, props & MEMORY), 1)
                ));
            else {
                context_ = ref_cnt_.find(file_key(filename, props & LARGE, props & MEMORY))->second.first;
                context_->grant(props & WRITE, props & REPLACE);
                ++ref_cnt_.find(file_key(filename, props & LARGE, props & MEMORY))->second.second;
            }
        }

        std::string archive::file_key(std::string filename, bool large, bool memory) const {
            return (large ? "l" : (memory ? "m" : "_")) + filename;
        }
    
#ifndef ALPS_NGS_SINGLE_THREAD
        boost::recursive_mutex archive::mutex_;
#endif
        std::map<std::string, std::pair<detail::archivecontext *, std::size_t> > archive::ref_cnt_;
    }
}
