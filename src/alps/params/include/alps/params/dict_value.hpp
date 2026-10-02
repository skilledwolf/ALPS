// Derived from ALPSCore params at 7146b9e1f017938a94e5dae35d88467cc5ba7969.
// Copyright (C) 1998-2018 ALPS Collaboration; modifications (C) 2026 ALPS Collaboration.
// SPDX-License-Identifier: MIT
// See ALPSCore-LICENSE.txt for the original permission notice.
#pragma once
#include <alps/params/dict_exceptions.hpp>
#include <alps/params_export.h>
#include <boost/serialization/complex.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/variant.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/variant.hpp>
#include <cmath>
#include <complex>
#include <cstdint>
#include <limits>
#include <ostream>
#include <string>
#include <type_traits>
#include <vector>

namespace alps::hdf5 {
class archive;
}
namespace alps::params_ns {
namespace detail {
struct None {
    template <class Archive> void serialize(Archive &, unsigned) {}
    bool operator==(const None &) const { return true; }
};
template <class T> struct vector_type : std::false_type {};
template <class T> struct vector_type<std::vector<T>> : std::true_type {
    using element = T;
};
template <class T>
struct allowed : std::bool_constant<std::is_arithmetic_v<T> || std::is_same_v<T, std::string> ||
                                    std::is_same_v<T, std::complex<double>>> {};
template <class T> struct allowed<std::vector<T>> : allowed<T> {};

// Core's checked signed/unsigned conversion, extended to reject lossy real
// conversion and to convert homogeneous vectors element by element.
template <class To, class From> To convert(const From &value, const std::string &key) {
    using namespace exception;
    if constexpr (std::is_same_v<To, From>)
        return value;
    else if constexpr (vector_type<To>::value && vector_type<From>::value) {
        To out;
        out.reserve(value.size());
        for (const auto &x : value)
            out.push_back(convert<typename vector_type<To>::element>(
                static_cast<typename vector_type<From>::element>(x), key));
        return out;
    } else if constexpr (std::is_integral_v<To> && std::is_integral_v<From> &&
                         !std::is_same_v<To, bool> && !std::is_same_v<From, bool>) {
        if constexpr (std::is_signed_v<From>) {
            if (value < 0) {
                if constexpr (!std::is_signed_v<To>)
                    throw value_mismatch(key, "negative value for unsigned integer");
                else if (value < std::numeric_limits<To>::lowest())
                    throw value_mismatch(key, "integer underflow");
                return static_cast<To>(value);
            }
        }
        if (static_cast<std::uintmax_t>(value) >
            static_cast<std::uintmax_t>(std::numeric_limits<To>::max()))
            throw value_mismatch(key, "integer overflow");
        return static_cast<To>(value);
    } else if constexpr (std::is_floating_point_v<To> && std::is_arithmetic_v<From> &&
                         !std::is_same_v<From, bool>) {
        const To out = static_cast<To>(value);
        // Compare in long double where available, and check integer precision
        // with binary digits so the rule also holds on Windows (double == long double).
        if constexpr (std::is_integral_v<From>) {
            if (!std::isfinite(out))
                throw value_mismatch(key, "real overflow");
            auto magnitude = static_cast<std::uintmax_t>(value);
            if constexpr (std::is_signed_v<From>)
                if (value < 0)
                    magnitude = std::uintmax_t(0) - magnitude;
            int bits = 0;
            for (auto n = magnitude; n; n >>= 1)
                ++bits;
            const int lost = bits - std::numeric_limits<To>::digits;
            if (lost > 0 && (magnitude & ((std::uintmax_t(1) << lost) - 1)))
                throw value_mismatch(key, "integer cannot be represented exactly as real");
        } else if (std::isfinite(value) && static_cast<From>(out) != value)
            throw value_mismatch(key, "lossy real conversion");
        return out;
    } else if constexpr (std::is_same_v<To, std::complex<double>> && std::is_arithmetic_v<From> &&
                         !std::is_same_v<From, bool>) {
        return {convert<double>(value, key), 0.};
    } else
        throw type_mismatch(key, "incompatible parameter types");
}

template <class T> auto canonical(const T &value, const std::string &key) {
    if constexpr (std::is_same_v<T, bool> || std::is_same_v<T, std::string> ||
                  std::is_same_v<T, std::complex<double>>)
        return value;
    else if constexpr (std::is_integral_v<T> && std::is_signed_v<T>)
        return convert<std::int64_t>(value, key);
    else if constexpr (std::is_integral_v<T>)
        return convert<std::uint64_t>(value, key);
    else if constexpr (std::is_floating_point_v<T>)
        return convert<double>(value, key);
    else if constexpr (vector_type<T>::value) {
        using E = typename vector_type<T>::element;
        using C = decltype(canonical(E{}, key));
        std::vector<C> out;
        out.reserve(value.size());
        for (const auto &x : value)
            out.push_back(canonical(static_cast<E>(x), key));
        return out;
    } else
        static_assert(allowed<T>::value, "Unsupported parameter type");
}
} // namespace detail

// Core's named value, empty state, assignment and as<T>/visitor API. Storage
// uses fixed-width families; parsing, interpreter objects and expressions do
// not belong to this class.
class ALPS_PARAMS_DECL dict_value {
  public:
    using None = detail::None;
    using value_type =
        boost::variant<None, bool, std::int64_t, std::uint64_t, double, std::complex<double>,
                       std::string, std::vector<bool>, std::vector<std::int64_t>,
                       std::vector<std::uint64_t>, std::vector<double>,
                       std::vector<std::complex<double>>, std::vector<std::string>>;
    dict_value() = default;
    explicit dict_value(const std::string &name) : name_(name) {}
    dict_value(const dict_value &) = default;
    dict_value &operator=(const dict_value &other) {
        val_ = other.val_;
        return *this;
    }
    bool empty() const { return val_.which() == 0; }
    void clear() { val_ = None{}; }
    template <class T> bool isType() const {
        return boost::apply_visitor(
            [](const auto &x) { return std::is_same_v<T, std::decay_t<decltype(x)>>; }, val_);
    }
    template <class T> const T &operator=(const T &value) {
        val_ = detail::canonical(value, name_);
        return value;
    }
    const char *operator=(const char *value) {
        val_ = std::string(value);
        return value;
    }
    template <class T> T as() const {
        static_assert(detail::allowed<T>::value, "Unsupported parameter type");
        if (empty())
            throw exception::uninitialized_value(name_, "value is not set");
        return boost::apply_visitor(
            [&](const auto &x) -> T { return detail::convert<T>(x, name_); }, val_);
    }
    template <class T, std::enable_if_t<detail::allowed<T>::value, int> = 0> operator T() const {
        return as<T>();
    }
    template <class F> decltype(auto) apply_visitor(const F &f) const {
        return boost::apply_visitor(f, val_);
    }
    bool equals(const dict_value &rhs) const { return val_ == rhs.val_; }
    void save(hdf5::archive &) const;
    void load(hdf5::archive &);
    friend ALPS_PARAMS_DECL std::ostream &operator<<(std::ostream &, const dict_value &);

  private:
    std::string name_;
    value_type val_;
    friend class boost::serialization::access;
    template <class Archive> void serialize(Archive &ar, unsigned) { ar & name_ & val_; }
};
template <class F> decltype(auto) apply_visitor(const F &f, const dict_value &v) {
    return v.apply_visitor(f);
}
} // namespace alps::params_ns
