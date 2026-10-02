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
} // namespace
void dict_value::save(hdf5::archive &ar) const {
    if (empty())
        throw exception::uninitialized_value(name_, "cannot checkpoint an unset value");
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
    // Decode the declared logical type, never guess from HDF5's physical type
    // (bool/int8 and complex shape are otherwise ambiguous).
    auto read = [&](auto exemplar) {
        using T = decltype(exemplar);
        if (type != logical_type<T>())
            return;
        auto validate = [&](auto element, bool array) {
            using E = decltype(element);
            constexpr bool complex = std::is_same_v<E, std::complex<double>>;
            using Physical = std::conditional_t<complex, double, E>;
            if (!ar.is_data("value") || !ar.is_datatype<Physical>("value") ||
                ar.is_complex("value") != complex)
                throw exception::type_mismatch(name_,
                                               "checkpoint payload disagrees with declared type");
            const auto shape = ar.extent("value");
            const bool valid =
                array ? (ar.is_null("value") ||
                         (shape.size() == (complex ? 2 : 1) && (!complex || shape.back() == 2)))
                      : (complex ? shape == std::vector<std::size_t>{2} : ar.is_scalar("value"));
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
        if (entry.second.empty())
            throw exception::uninitialized_value(entry.first, "cannot checkpoint an unset value");
    // Index entries so parameter names containing '/' or TOML punctuation do
    // not become archive paths. Overwriting also removes stale prior entries.
    if (ar.is_group("entries"))
        ar.delete_group("entries");
    ar.create_group("entries");
    ar["format"] << std::string("alps.params.v1");
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
    if (format != "alps.params.v1")
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
