// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/complex.hpp>
#include <alps/hdf5/vector.hpp>
#include <alps/params.hpp>
namespace alps::params_ns {
namespace {
void validate_checkpoint(const dict_value &value, const std::string &key) {
    // The archive stores C strings. Reject unrepresentable names and values
    // before overwriting a checkpoint rather than silently truncating them.
    if (key.find('\0') != std::string::npos)
        throw std::invalid_argument("Cannot checkpoint a parameter name containing NUL");
    if (value.empty())
        throw exception::uninitialized_value(key, "cannot checkpoint an unset value");
    auto validate_string = [&](const std::string &text) {
        if (text.find('\0') != std::string::npos)
            throw exception::value_mismatch(key, "cannot checkpoint a string containing NUL");
    };
    value.apply_visitor([&](const auto &text) {
        using T = std::decay_t<decltype(text)>;
        if constexpr (std::is_same_v<T, std::string>)
            validate_string(text);
        else if constexpr (std::is_same_v<T, std::vector<std::string>>)
            for (const auto &element : text)
                validate_string(element);
    });
}
} // namespace
void dict_value::save(hdf5::archive &ar) const {
    validate_checkpoint(*this, name_);
    apply_visitor([&](const auto &value) {
        using T = std::decay_t<decltype(value)>;
        if constexpr (!std::is_same_v<T, None>)
            ar[ar.get_context()] << value;
    });
}
void dict_value::load(hdf5::archive &ar) {
    const auto path = ar.get_context();
    if (!ar.is_data(path) || ar.is_null(path))
        throw exception::type_mismatch(name_, "checkpoint value must be a native dataset");
    const auto shape = ar.extent(path);
    if (shape.size() > 1)
        throw exception::type_mismatch(name_, "checkpoint values must be scalars or vectors");
    value_type loaded;
    bool matched = false;
    // Canonical HDF5 datatype and rank identify the variant without a second
    // type tag that could disagree with the payload.
    auto read = [&](auto exemplar) {
        using T = decltype(exemplar);
        auto matches = [&](auto element, bool array) {
            using E = decltype(element);
            return shape.size() == (array ? 1 : 0) && ar.is_datatype<E>(path);
        };
        bool compatible;
        if constexpr (detail::vector_type<T>::value)
            compatible = matches(typename T::value_type{}, true);
        else
            compatible = matches(T{}, false);
        if (!compatible)
            return;
        T v;
        ar[path] >> v;
        loaded = std::move(v);
        matched = true;
    };
    read(bool{});
    read(std::int64_t{});
    read(std::uint64_t{});
    read(double{});
    read(std::complex<double>{});
    read(std::string{});
    read(std::vector<bool>{});
    read(std::vector<std::int64_t>{});
    read(std::vector<std::uint64_t>{});
    read(std::vector<double>{});
    read(std::vector<std::complex<double>>{});
    read(std::vector<std::string>{});
    if (!matched)
        throw exception::type_mismatch(name_, "unsupported checkpoint datatype");
    val_.swap(loaded);
}
void dictionary::save(hdf5::archive &ar) const {
    for (const auto &entry : *this)
        validate_checkpoint(entry.second, entry.first);
    // Index entries so parameter names containing '/' or TOML punctuation do
    // not become archive paths. Overwriting also removes stale prior entries.
    if (ar.is_group("entries"))
        ar.delete_group("entries");
    ar.create_group("entries");
    ar["format"] << std::string("alps.params.v2");
    std::size_t i = 0;
    for (const auto &entry : *this) {
        const auto path = "entries/" + std::to_string(i++);
        ar[path + "/name"] << entry.first;
        ar[path + "/value"] << entry.second;
    }
}
void dictionary::load(hdf5::archive &ar) {
    std::string format;
    if (ar.is_data("format"))
        ar["format"] >> format;
    if (format != "alps.params.v2")
        throw std::runtime_error(
            "Unsupported params checkpoint format; legacy checkpoints require offline conversion");
    dictionary loaded;
    const auto entries = ar.list_children("entries");
    for (std::size_t i = 0; i < entries.size(); ++i) {
        const auto path = "entries/" + std::to_string(i);
        std::string name;
        ar[path + "/name"] >> name;
        if (loaded.exists(name))
            throw std::runtime_error("Duplicate parameter in checkpoint: " + name);
        ar[path + "/value"] >> loaded[name];
    }
    swap(*this, loaded);
}
} // namespace alps::params_ns
