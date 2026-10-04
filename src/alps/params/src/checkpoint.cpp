// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/complex.hpp>
#include <alps/hdf5/vector.hpp>
#include <alps/params.hpp>
namespace alps::params_ns {
namespace {
template <class T> const char *type_name();
#define ALPS_PARAM_TYPE(T, NAME)                                                                   \
    template <> const char *type_name<T>() { return NAME; }
ALPS_PARAM_TYPE(bool, "bool")
ALPS_PARAM_TYPE(std::int64_t, "int64")
ALPS_PARAM_TYPE(std::uint64_t, "uint64")
ALPS_PARAM_TYPE(double, "float64")
ALPS_PARAM_TYPE(std::complex<double>, "complex128")
ALPS_PARAM_TYPE(std::string, "string")
#undef ALPS_PARAM_TYPE
template <class T> std::string logical_type() {
    if constexpr (detail::vector_type<T>::value)
        return std::string(type_name<typename detail::vector_type<T>::element>()) + "[]";
    else
        return type_name<T>();
}
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
        if constexpr (!std::is_same_v<T, None>) {
            ar["type"] << logical_type<T>();
            ar["value"] << value;
        }
    });
}
void dict_value::load(hdf5::archive &ar) {
    std::string type;
    ar["type"] >> type;
    value_type loaded;
    bool matched = false;
    // The logical type and native HDF5 datatype must agree exactly.
    auto read = [&](auto exemplar) {
        using T = decltype(exemplar);
        if (type != logical_type<T>())
            return;
        auto validate = [&](auto element, bool array) {
            using E = decltype(element);
            if (!ar.is_data("value") || !ar.is_datatype<E>("value"))
                throw exception::type_mismatch(name_,
                                               "checkpoint payload disagrees with declared type");
            const auto shape = ar.extent("value");
            const bool valid = !ar.is_null("value") &&
                               (array ? shape.size() == 1 : ar.is_scalar("value"));
            if (!valid)
                throw exception::type_mismatch(name_, "invalid checkpoint shape");
        };
        if constexpr (detail::vector_type<T>::value)
            validate(typename T::value_type{}, true);
        else
            validate(T{}, false);
        T v;
        ar["value"] >> v;
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
        throw exception::type_mismatch(name_, "unknown checkpoint type '" + type + "'");
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
        ar[path] << entry.second;
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
        ar[path] >> loaded[name];
    }
    swap(*this, loaded);
}
} // namespace alps::params_ns
