// Copyright (C) 1998-2018, 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#pragma once
#include <alps/alea/core.hpp>
#include <alps/hdf5/archive.hpp>
#include <algorithm>

namespace alps::alea {
// One ndview boundary for native canonical HDF5 primitives. Entering groups
// changes only the logical path; reading never creates or modifies objects.
class hdf5_serializer : public serializer, public deserializer {
public:
    hdf5_serializer(hdf5::archive& archive, std::string path)
        : archive_(archive), path_(archive.complete_path(std::move(path))) {}
    void enter(std::string const& group) override {
        if (!group.empty()) validate_key(group); // Empty group denotes the current path.
        groups_.push_back(group);
    }
    void exit() override {
        if (groups_.empty()) throw std::runtime_error("exit without enter");
        groups_.pop_back();
    }
    std::vector<size_t> get_shape(std::string const& key) override {
        auto path = get_path(key);
        if (archive_.is_null(path)) throw hdf5::wrong_type("NULL has no ALEA array shape");
        return archive_.extent(path);
    }
#define ALPS_ALEA_PRIMITIVE(T) \
    void write(std::string const& key, ndview<T const> data) override { do_write(key, data); } \
    void read(std::string const& key, ndview<T> data) override { do_read(key, data); }
    ALPS_ALEA_PRIMITIVE(double)
    ALPS_ALEA_PRIMITIVE(std::complex<double>)
    ALPS_ALEA_PRIMITIVE(int64_t)
    ALPS_ALEA_PRIMITIVE(uint64_t)
    ALPS_ALEA_PRIMITIVE(int32_t)
    ALPS_ALEA_PRIMITIVE(uint32_t)
#undef ALPS_ALEA_PRIMITIVE
private:
    static void validate_key(std::string const& key) {
        if (key.empty() || key == "." || key == ".." || key.find('/') != std::string::npos)
            throw std::runtime_error("ALEA serialization keys must be single nonempty segments");
    }
    std::string get_path(std::string const& key) const {
        std::string result = path_;
        for (auto const& group : groups_) result += "/" + group;
        if (!key.empty()) { validate_key(key); result += "/" + key; }
        return result;
    }
    template<class T> void do_write(std::string const& key, ndview<T const> data) {
        auto const group = get_path("");
        if (!archive_.is_group(group)) archive_.create_group(group);
        auto const path = get_path(key);
        if (!data.ndim()) {
            if (!data.data()) throw std::runtime_error("missing ALEA scalar buffer");
            archive_.write(path, *data.data());
        } else {
            std::vector<size_t> shape(data.shape(), data.shape() + data.ndim());
            archive_.write(path, data.data(), shape);
        }
    }
    template<class T> void do_read(std::string const& key, ndview<T> data) {
        auto shape = get_shape(key);
        if (shape.size() != data.ndim() || !std::equal(shape.begin(), shape.end(), data.shape()))
            throw size_mismatch();
        if (!archive_.is_datatype<T>(get_path(key)))
            throw hdf5::wrong_type("unexpected ALEA primitive datatype at " + get_path(key));
        // A null destination deliberately discards a derived field after
        // validating its shape and datatype. Empty arrays have no payload to transfer.
        if (!data.data() || !data.size()) return;
        auto const path = get_path(key);
        if (shape.empty()) archive_.read(path, *data.data());
        else archive_.read(path, data.data(), shape);
    }
    hdf5::archive& archive_;
    std::string path_;
    std::vector<std::string> groups_;
};
}
