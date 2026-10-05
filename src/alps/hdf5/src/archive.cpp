// Copyright (C) 2010-2012 Lukas Gamper; ALPS contributors.
// SPDX-License-Identifier: MIT
#include <alps/hdf5/archive.hpp>

#include <highfive/highfive.hpp>
#include <boost/filesystem.hpp>
#include <algorithm>
#include <exception>
#include <iostream>
#include <limits>
#include <mutex>
#include <type_traits>

namespace alps::hdf5 {
namespace detail {
namespace hf = HighFive;

// Serialize native operations for non-threadsafe HDF5 providers. A caller must
// still serialize save/load transactions using the same archive's context.
std::recursive_mutex mutex;

struct native_type : hf::DataType {
    explicit native_type(hf::DataType const& stored)
        : hf::DataType(H5Tget_native_type(stored.getId(), H5T_DIR_ASCEND)) {
        if (!isValid()) throw wrong_type("cannot inspect HDF5 datatype");
    }
};

template<class T> bool matches(hf::DataType const& stored) {
    if constexpr (std::is_same_v<T, std::string>)
        return stored.getClass() == hf::DataTypeClass::String;
    else
        return native_type(stored) == hf::create_datatype<T>();
}

bool complex_type(hf::DataType const& type) {
    return type.getClass() == hf::DataTypeClass::Compound &&
        (matches<std::complex<float>>(type) || matches<std::complex<double>>(type) ||
         matches<std::complex<long double>>(type));
}

struct file_handle : hf::File {
    using hf::File::File;

    // HighFive's destructor is a nonthrowing fallback. Explicit finalization
    // must report close failures before a caller publishes a checkpoint.
    void close() {
        if (H5Fclose(_hid) < 0) throw archive_error("cannot close HDF5 file " + getName());
        _hid = H5I_INVALID_HID;
    }
};

struct archivecontext {
    std::string const filename;
    bool const writable;
    file_handle file;

    archivecontext(std::string name, std::string const& mode)
        : filename(std::move(name)), writable(mode != "r"),
          file(filename, mode == "r" ? hf::File::ReadOnly :
                         mode == "a" ? hf::File::OpenOrCreate : hf::File::Truncate) {}
};

// Called under mutex. Finalize every view of this file, even if flushing fails;
// HighFive RAII remains the fallback if checked close itself cannot succeed.
void finish(std::shared_ptr<archivecontext> const& context, bool flush) {
    if (!context || !context->file.isValid()) return;
    hf::SilenceHDF5 silence;
    std::exception_ptr failure;
    if (flush) {
        try { context->file.flush(); }
        catch (...) { failure = std::current_exception(); }
    }
    try { context->file.close(); }
    catch (...) { if (!failure) failure = std::current_exception(); }
    if (failure) {
        try { std::rethrow_exception(failure); }
        catch (hf::Exception const& error) { throw archive_error(error.what()); }
    }
}

template<class F> decltype(auto) access(std::shared_ptr<archivecontext> const& context, F&& fn) {
    std::lock_guard<std::recursive_mutex> lock(mutex);
    if (!context || !context->file.isValid()) throw archive_closed("the archive is closed");
    hf::SilenceHDF5 silence;
    try { return fn(*context); }
    catch (hf::Exception const& error) { throw archive_error(error.what()); }
}

struct location {
    std::string object, attribute;
    bool is_attribute;
    explicit location(std::string const& path) : object(path), is_attribute(false) {
        auto at = path.find_last_of('@');
        if (at != std::string::npos) {
            is_attribute = true;
            attribute = path.substr(at + 1);
            object = path.substr(0, at);
            while (object.size() > 1 && object.back() == '/') object.pop_back();
            if (object.empty()) object = "/";
            if (attribute.empty()) throw invalid_path("empty attribute name: " + path);
        }
    }
};

template<class F> decltype(auto) node(archivecontext& context, std::string const& path, F&& fn) {
    auto& file = context.file;
    if (!file.exist(path)) throw path_not_found("no object at " + path);
    if (file.getObjectType(path) == hf::ObjectType::Group) {
        auto group = file.getGroup(path);
        return fn(group);
    }
    if (file.getObjectType(path) == hf::ObjectType::Dataset) {
        auto dataset = file.getDataSet(path);
        return fn(dataset);
    }
    throw wrong_type("expected group or dataset at " + path);
}

template<class F> decltype(auto) stored(archivecontext& context, location const& path, F&& fn) {
    if (path.is_attribute)
        return node(context, path.object, [&](auto& parent) -> decltype(auto) {
            if (!parent.hasAttribute(path.attribute))
                throw path_not_found("no attribute " + path.attribute + " on " + path.object);
            auto attr = parent.getAttribute(path.attribute);
            return fn(attr);
        });
    if (!context.file.exist(path.object)) throw path_not_found("no dataset at " + path.object);
    if (context.file.getObjectType(path.object) != hf::ObjectType::Dataset)
        throw wrong_type("expected dataset at " + path.object);
    auto dataset = context.file.getDataSet(path.object);
    return fn(dataset);
}

size_t product(std::vector<size_t> const& dimensions, size_t value = 1) {
    if (std::find(dimensions.begin(), dimensions.end(), 0) != dimensions.end()) return 0;
    for (auto n : dimensions) {
        if (n > std::numeric_limits<size_t>::max() / value) throw archive_error("HDF5 extent product overflows");
        value *= n;
    }
    return value;
}

void selection(std::vector<size_t> const& size, std::vector<size_t> const& count,
               std::vector<size_t> const& offset) {
    if (size.size() != count.size() || size.size() != offset.size())
        throw archive_error("selection rank does not match dataset rank");
    for (size_t i = 0; i < size.size(); ++i)
        if (offset[i] > size[i] || count[i] > size[i] - offset[i])
            throw archive_error("selection exceeds dataset extent");
}

template<class T> void check_read_type(hf::DataType const& type) {
    using C = hf::DataTypeClass;
    auto expected = hf::create_datatype<T>().getClass();
    auto actual = type.getClass();
    bool compatible = expected == C::String ? actual == C::String :
        expected == C::Enum ? matches<bool>(type) :
        expected == C::Compound ? complex_type(type) : actual == C::Integer || actual == C::Float;
    if (!compatible) throw wrong_type("incompatible HDF5 datatype; convert legacy archives with alps-hdf5-convert");
}

// HighFive's generic string reader reclaims after unserializing. Own that memory
// through string assignment too, so allocation failures cannot leak HDF5 buffers.
struct string_memory {
    hf::DataType type;
    hf::DataSpace space;
    std::vector<char*>& pointers;
    ~string_memory() { H5Treclaim(type.getId(), space.getId(), H5P_DEFAULT, pointers.data()); }
};

template<class T, class Object>
void read_values(Object const& object, hf::DataType const& type, T* values, size_t count) {
    if constexpr (std::is_same_v<T, std::string>) {
        if (type.isVariableStr()) {
            std::vector<char*> raw(count, nullptr);
            string_memory cleanup{type, hf::DataSpace(std::vector<size_t>{count}), raw};
            object.read_raw(raw.data(), type);
            for (size_t i = 0; i < count; ++i) values[i] = raw[i] ? raw[i] : "";
        } else {
            auto width = type.getSize();
            product({count}, width);
            std::vector<char> raw(count * width);
            object.read_raw(raw.data(), type);
            auto padding = type.asStringType().getPadding();
            for (size_t i = 0; i < count; ++i) {
                auto begin = raw.data() + i * width, end = begin + width;
                if (padding == hf::StringPadding::SpacePadded)
                    while (end != begin && end[-1] == ' ') --end;
                else end = std::find(begin, end, '\0');
                values[i].assign(begin, end);
            }
        }
    } else if constexpr (std::is_same_v<T, bool>) {
        // Read into bytes first: malformed enum values must never become an
        // invalid C++ bool representation, even transiently.
        std::vector<unsigned char> raw(count);
        static_assert(sizeof(bool) == sizeof(unsigned char));
        object.read_raw(raw.data(), hf::create_datatype<bool>());
        if (std::any_of(raw.begin(), raw.end(), [](auto value) { return value > 1; }))
            throw wrong_type("Boolean dataset contains a value other than FALSE or TRUE");
        for (size_t i = 0; i < count; ++i) values[i] = raw[i] != 0;
    } else {
        object.read_raw(values, hf::create_datatype<T>());
    }
}

template<class T>
void read(archivecontext& context, location const& path, T* values, bool scalar,
          std::vector<size_t> count = {}, std::vector<size_t> offset = {}) {
    stored(context, path, [&](auto& object) {
        auto space = object.getSpace();
        auto size = space.getDimensions();
        auto kind = H5Sget_simple_extent_type(space.getId());
        if (kind == H5S_NULL) throw wrong_type("NULL dataspace has no array shape; convert it using its application schema");
        if ((kind == H5S_SCALAR) != scalar) throw wrong_type("scalar/array shape mismatch");
        auto type = object.getDataType();
        check_read_type<T>(type);
        if (!scalar) {
            if (offset.empty()) offset.resize(size.size(), 0);
            selection(size, count, offset);
            if (path.is_attribute && count != size) throw std::logic_error("partial attribute reads are unsupported");
        }
        auto n = scalar ? 1 : product(count, sizeof(T)) / sizeof(T);
        if (!n) return;
        if (!values) throw archive_error("null read buffer for a nonempty selection");
        if constexpr (std::is_same_v<std::decay_t<decltype(object)>, hf::DataSet>) {
            if (!scalar) { read_values(object.select(offset, count), type, values, n); return; }
        }
        read_values(object, type, values, n);
    });
}

template<class T>
void write(archivecontext& context, location const& path, T const* values, bool scalar,
           std::vector<size_t> size = {}, std::vector<size_t> count = {}, std::vector<size_t> offset = {}) {
    if (!context.writable) throw archive_error("the archive is not writable");
    if (!scalar && size.empty()) throw archive_error("array writes require an explicit shape, including empty arrays");
    if (!scalar) {
        if (count.empty()) count = size;
        if (offset.empty()) offset.resize(size.size(), 0);
        selection(size, count, offset);
        product(size, sizeof(T));
        if (path.is_attribute && count != size) throw std::logic_error("partial attribute writes are unsupported");
    }
    auto n = scalar ? 1 : product(count, sizeof(T)) / sizeof(T);
    if (n && !values) throw archive_error("null write buffer for a nonempty selection");
    std::vector<char const*> strings;
    if constexpr (std::is_same_v<T, std::string>) {
        strings.reserve(n);
        for (size_t i = 0; i < n; ++i) {
            if (values[i].find('\0') != std::string::npos) throw archive_error("HDF5 strings cannot contain NUL bytes");
            strings.push_back(values[i].c_str());
        }
    }
    auto type = hf::create_datatype<T>();
    auto space = scalar ? hf::DataSpace::Scalar() : hf::DataSpace(size);
    auto matches_target = [&](auto const& object) {
        // H5Tequal ignores the string character set, but HDF5 cannot convert
        // ASCII storage to UTF-8 during a write. Replace that storage explicitly.
        if constexpr (std::is_same_v<T, std::string>) {
            auto stored = object.getDataType();
            if (stored.getClass() != hf::DataTypeClass::String ||
                stored.asStringType().getCharacterSet() != type.asStringType().getCharacterSet())
                return false;
        }
        return object.getDataType() == type && object.getSpace().getDimensions() == size &&
            H5Sget_simple_extent_type(object.getSpace().getId()) == (scalar ? H5S_SCALAR : H5S_SIMPLE);
    };
    auto transfer = [&](auto& target) {
        if (!n) return;
        try {
            if constexpr (std::is_same_v<T, std::string>) target.write_raw(strings.data(), type);
            else target.write_raw(values, type);
        } catch (hf::Exception const& error) {
            throw archive_error(path.object + (path.is_attribute ? "/@" + path.attribute : "") + ": " + error.what());
        }
    };
    if (path.is_attribute) {
        node(context, path.object, [&](auto& parent) {
            if (parent.hasAttribute(path.attribute) && !matches_target(parent.getAttribute(path.attribute)))
                parent.deleteAttribute(path.attribute);
            auto attr = parent.hasAttribute(path.attribute) ? parent.getAttribute(path.attribute) :
                parent.createAttribute(path.attribute, space, type);
            transfer(attr);
        });
    } else {
        auto& file = context.file;
        if (file.exist(path.object)) {
            bool same = file.getObjectType(path.object) == hf::ObjectType::Dataset &&
                matches_target(file.getDataSet(path.object));
            if (!same) file.unlink(path.object);
        }
        auto dataset = file.exist(path.object) ? file.getDataSet(path.object) :
            file.createDataSet(path.object, space, type);
        if (scalar) transfer(dataset);
        else if (n) { auto selected = dataset.select(offset, count); transfer(selected); }
    }
}
} // namespace detail

archive::archive(boost::filesystem::path const& filename, std::string mode) {
    std::lock_guard<std::recursive_mutex> lock(detail::mutex);
    if (mode != "r" && mode != "a" && mode != "w")
        throw archive_error("invalid HDF5 open mode: " + mode);
    auto name = boost::filesystem::absolute(filename).lexically_normal().string();
    detail::hf::SilenceHDF5 silence;
    try {
        context_ = std::make_shared<detail::archivecontext>(name, mode);
    } catch (detail::hf::Exception const& error) {
        throw archive_not_found(error.what());
    }
}
archive::archive(archive const& other) : current_(other.current_), context_(other.context_) {}
archive::~archive() {
    std::lock_guard<std::recursive_mutex> lock(detail::mutex);
    if (context_.use_count() == 1) {
        try { detail::finish(context_, true); }
        catch (std::exception const& error) { std::cerr << "Closing archive: " << error.what() << '\n'; }
    }
    context_.reset();
}
void archive::close() {
    std::lock_guard<std::recursive_mutex> lock(detail::mutex);
    detail::finish(context_, true);
    context_.reset();
}
bool archive::is_open() {
    std::lock_guard<std::recursive_mutex> lock(detail::mutex);
    return context_ && context_->file.isValid();
}

void save_checkpoint(boost::filesystem::path const& filename,
                     std::function<void(archive&)> const& save) {
    archive::publish(filename,save,false);
}
void update_archive(boost::filesystem::path const& filename,
                    std::function<void(archive&)> const& update) {
    archive::publish(filename,update,true);
}
void archive::publish(boost::filesystem::path const& filename,
                      std::function<void(archive&)> const& save,bool preserve) {
    auto const target = boost::filesystem::absolute(filename).lexically_normal();
    boost::filesystem::path temporary_directory;
    // Reserving a private directory avoids truncating an existing temporary
    // file, even if another save chooses the same random name concurrently.
    do {
        temporary_directory = target.parent_path() / boost::filesystem::unique_path(
            target.filename().string() + ".tmp.%%%%-%%%%-%%%%-%%%%");
    } while (!boost::filesystem::create_directory(temporary_directory));
    try {
        auto const temporary = temporary_directory / "checkpoint.h5";
        {
            if (preserve) boost::filesystem::copy_file(target,temporary);
            archive ar(temporary, preserve ? "a" : "w");
            // Retain finalization ownership even if the callback closes its
            // handle or keeps a copy. Escaped handles cannot defer publication.
            auto context = ar.context_;
            std::exception_ptr failure;
            try { save(ar); }
            catch (...) { failure = std::current_exception(); }
            {
                std::lock_guard<std::recursive_mutex> lock(detail::mutex);
                try { detail::finish(context, !failure); }
                catch (...) { if (!failure) failure = std::current_exception(); }
                ar.context_.reset();
                context.reset();
            }
            if (failure) std::rethrow_exception(failure);
        }
        boost::filesystem::rename(temporary, target);
    } catch (...) {
        boost::system::error_code cleanup_error;
        boost::filesystem::remove_all(temporary_directory, cleanup_error);
        throw;
    }
    boost::system::error_code cleanup_error;
    boost::filesystem::remove(temporary_directory, cleanup_error);
}

std::string const& archive::get_filename() const {
    if (!context_) throw archive_closed("the archive is closed");
    return context_->filename;
}
std::string archive::get_context() const { return current_; }
void archive::set_context(std::string const& path) { current_ = complete_path(path); }
std::string archive::complete_path(std::string path) const {
    if (path.empty() || path.front() != '/') path = current_ + (current_ == "/" || path.empty() ? "" : "/") + path;
    std::vector<std::string> segments;
    for (size_t start = 0, end; start <= path.size(); start = end + 1) {
        end = path.find('/', start);
        if (end == std::string::npos) end = path.size();
        auto part = path.substr(start, end - start);
        if (part == "..") { if (!segments.empty()) segments.pop_back(); }
        else if (!part.empty() && part != ".") segments.push_back(std::move(part));
    }
    std::string result;
    for (auto const& part : segments) result += "/" + part;
    return result.empty() ? "/" : result;
}
std::string archive::encode_segment(std::string text) const {
    for (auto c : {'&', '/'}) {
        size_t pos = 0;
        auto replacement = "&#" + std::to_string(int(c)) + ";";
        while ((pos = text.find(c, pos)) != std::string::npos) { text.replace(pos, 1, replacement); pos += replacement.size(); }
    }
    return text;
}
std::string archive::decode_segment(std::string text) const {
    size_t pos = 0;
    while ((pos = text.find("&#", pos)) != std::string::npos) {
        if (text.compare(pos, 5, "&#47;") == 0) text.replace(pos, 5, "/");
        else if (text.compare(pos, 5, "&#38;") == 0) text.replace(pos, 5, "&");
        ++pos;
    }
    return text;
}

bool archive::is_data(std::string path) const {
    detail::location target(complete_path(path));
    if (target.is_attribute) throw invalid_path("not a dataset path: " + path);
    return detail::access(context_, [&](auto& c) { return c.file.exist(target.object) && c.file.getObjectType(target.object) == detail::hf::ObjectType::Dataset; });
}
bool archive::is_group(std::string path) const {
    detail::location target(complete_path(path));
    return detail::access(context_, [&](auto& c) { return !target.is_attribute && c.file.exist(target.object) && c.file.getObjectType(target.object) == detail::hf::ObjectType::Group; });
}
bool archive::is_attribute(std::string path) const {
    detail::location target(complete_path(path));
    return detail::access(context_, [&](auto& c) {
        return target.is_attribute && c.file.exist(target.object) && detail::node(c, target.object, [&](auto& object) { return object.hasAttribute(target.attribute); });
    });
}
bool archive::is_scalar(std::string path) const {
    return detail::access(context_, [&](auto& c) { return detail::stored(c, detail::location(complete_path(path)), [](auto& object) { return H5Sget_simple_extent_type(object.getSpace().getId()) == H5S_SCALAR; }); });
}
bool archive::is_null(std::string path) const {
    return detail::access(context_, [&](auto& c) { return detail::stored(c, detail::location(complete_path(path)), [](auto& object) { return H5Sget_simple_extent_type(object.getSpace().getId()) == H5S_NULL; }); });
}
bool archive::is_complex(std::string path) const {
    if (is_group(path)) return false;
    return detail::access(context_, [&](auto& c) { return detail::stored(c, detail::location(complete_path(path)), [](auto& object) { return detail::complex_type(object.getDataType()); }); });
}
std::vector<size_t> archive::extent(std::string path) const {
    return detail::access(context_, [&](auto& c) { return detail::stored(c, detail::location(complete_path(path)), [](auto& object) { return object.getSpace().getDimensions(); }); });
}
size_t archive::dimensions(std::string path) const { return extent(path).size(); }
std::vector<std::string> archive::list_children(std::string path) const {
    return detail::access(context_, [&](auto& c) { return c.file.getGroup(complete_path(path)).listObjectNames(); });
}
std::vector<std::string> archive::list_attributes(std::string path) const {
    detail::location target(complete_path(path));
    if (target.is_attribute) throw invalid_path("attributes cannot have attributes");
    return detail::access(context_, [&](auto& c) { return detail::node(c, target.object, [](auto& object) { return object.listAttributeNames(); }); });
}
void archive::create_group(std::string path) const {
    detail::location target(complete_path(path));
    if (target.is_attribute) throw invalid_path("not a group path: " + path);
    detail::access(context_, [&](auto& c) {
        if (!c.writable) throw archive_error("the archive is not writable");
        if (c.file.exist(target.object)) {
            if (c.file.getObjectType(target.object) == detail::hf::ObjectType::Group) return;
            c.file.unlink(target.object);
        }
        c.file.createGroup(target.object);
    });
}
void archive::delete_data(std::string path) const {
    if (is_group(path)) throw invalid_path("dataset path contains a group");
    detail::access(context_, [&](auto& c) { if (!c.writable) throw archive_error("the archive is not writable"); if (is_data(path)) c.file.unlink(complete_path(path)); });
}
void archive::delete_group(std::string path) const {
    if (is_data(path)) throw invalid_path("group path contains a dataset");
    detail::access(context_, [&](auto& c) { if (!c.writable) throw archive_error("the archive is not writable"); if (is_group(path)) c.file.unlink(complete_path(path)); });
}
void archive::delete_attribute(std::string path) const {
    detail::location target(complete_path(path));
    if (!target.is_attribute) throw invalid_path("not an attribute path: " + path);
    detail::access(context_, [&](auto& c) {
        if (!c.writable) throw archive_error("the archive is not writable");
        detail::node(c, target.object, [&](auto& parent) { if (parent.hasAttribute(target.attribute)) parent.deleteAttribute(target.attribute); });
    });
}
detail::archive_proxy<archive> archive::operator[](std::string const& path) { return detail::archive_proxy<archive>(path, *this); }

#define ALPS_DEFINE_IO(T) \
void archive::read(std::string path, T& value) const { detail::access(context_, [&](auto& c) { detail::read(c, detail::location(complete_path(path)), &value, true); }); } \
void archive::read(std::string path, T* value, std::vector<size_t> count, std::vector<size_t> offset) const { detail::access(context_, [&](auto& c) { detail::read(c, detail::location(complete_path(path)), value, false, count, offset); }); } \
void archive::write(std::string path, T value) const { detail::access(context_, [&](auto& c) { detail::write(c, detail::location(complete_path(path)), &value, true); }); } \
void archive::write(std::string path, T const* value, std::vector<size_t> size, std::vector<size_t> count, std::vector<size_t> offset) const { detail::access(context_, [&](auto& c) { detail::write(c, detail::location(complete_path(path)), value, false, size, count, offset); }); } \
bool archive::is_datatype_impl(std::string path, T) const { return detail::access(context_, [&](auto& c) { return detail::stored(c, detail::location(complete_path(path)), [](auto& object) { return detail::matches<T>(object.getDataType()); }); }); } \
namespace detail { \
T* get_pointer<T>::apply(T& value) { return &value; } \
T const* get_pointer<T const>::apply(T const& value) { return &value; } \
bool is_vectorizable<T>::apply(T const&) { return true; } \
bool is_vectorizable<T const>::apply(T&) { return true; } \
} \
void save(archive& ar, std::string const& path, T const& value, std::vector<size_t> size, std::vector<size_t> count, std::vector<size_t> offset) { if (size.empty()) ar.write(path, value); else ar.write(path, &value, size, count, offset); } \
void load(archive& ar, std::string const& path, T& value, std::vector<size_t> count, std::vector<size_t> offset) { if (count.empty()) ar.read(path, value); else ar.read(path, &value, count, offset); }
ALPS_NGS_FOREACH_NATIVE_HDF5_TYPE(ALPS_DEFINE_IO)
#undef ALPS_DEFINE_IO
} // namespace alps::hdf5
