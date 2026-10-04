/*
 * Copyright (C) 1998-2018 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 * For use in publications, see ACKNOWLEDGE.TXT
 */
#pragma once

#include <alps/alea/core.hpp>
#include <alps/serialization/core.hpp>
#include <alps/serialization/util.hpp>
#include <alps/serialization/eigen.hpp>

#include <array>

namespace alps { namespace serialization {

/** Serializes Eigen array of complex_op<double> */
template <typename Derived>
typename std::enable_if<
        eigen_scalar_is<Derived, alps::alea::complex_op<double>>::value>::type
serialize(serializer &ser, const std::string &key,
          const Eigen::PlainObjectBase<Derived> &value)
{
    using scalar_type = eigen_scalar_t<Derived>;
    using matrix_type = Eigen::Matrix<scalar_type, Eigen::Dynamic, Eigen::Dynamic>;

    // Ensure that evaluated expression will be Fortran-contiguous
    if (!eigen_is_contiguous(value))
        return serialize(ser, key, matrix_type(value));

    // complex_op<...> is a wrapper around double[2][2].  In order to convert
    // it to serializable object, we first cast it to doubles, and then add
    // two further dimensions.
    const double *dbl_data = reinterpret_cast<const double *>(value.data());
    if (Derived::ColsAtCompileTime == 1 || Derived::RowsAtCompileTime == 1) {
        // Omit second dimension for simple vectors
        std::array<size_t, 3> dims = {{(size_t)value.size(), 2, 2}};
        ser.write(key, ndview<const double>(dbl_data, dims.data(), dims.size()));
    } else {
        // Eigen arrays are column-major
        std::array<size_t, 4> dims = {{(size_t)value.cols(), (size_t)value.rows(), 2, 2}};
        ser.write(key, ndview<const double>(dbl_data, dims.data(), dims.size()));
    }
}

/** Deserializes Eigen array or matrix of complex_op<double> */
template <typename Derived>
typename std::enable_if<
        eigen_scalar_is<Derived, alps::alea::complex_op<double>>::value>::type
deserialize(deserializer &ser, const std::string &key,
            Eigen::PlainObjectBase<Derived> &value)
{
    // Ensure that evaluated expression will be Fortran-contiguous
    if (!eigen_is_contiguous(value)) {
        using matrix_type = Eigen::Matrix<eigen_scalar_t<Derived>, Eigen::Dynamic, Eigen::Dynamic>;
        matrix_type staged(value.rows(), value.cols());
        deserialize(ser, key, staged);
        value = staged;
        return;
    }

    // complex_op<...> is a wrapper around double[2][2].  In order to convert
    // it to serializable object, we first cast it to doubles, and then add
    // two further dimensions.
    double *dbl_data = reinterpret_cast<double *>(value.data());
    if (Derived::ColsAtCompileTime == 1 || Derived::RowsAtCompileTime == 1) {
        // Omit second dimension for simple vectors
        std::array<size_t, 3> shape = {{(size_t)value.size(), 2, 2}};
        ser.read(key, ndview<double>(dbl_data, shape.data(), shape.size()));
    } else {
        // Extract underlying buffer and read
        std::array<size_t, 4> shape = {{(size_t)value.cols(), (size_t)value.rows(), 2, 2}};
        ser.read(key, ndview<double>(dbl_data, shape.data(), shape.size()));
    }
}

}}

namespace alps::alea::internal {
enum class result_kind : uint32_t { mean=1, variance=2, covariance=3, autocorrelation=4, batch=5, batch_accumulator=6 };
struct result_sentry : serializer_sentry {
    result_sentry(serializer& s, std::string const& key, result_kind kind)
        : serializer_sentry(s, key) {
        alps::serialization::serialize(s, "@version", uint64_t(1));
        alps::serialization::serialize(s, "@kind", static_cast<uint32_t>(kind));
    }
};
struct result_reader_sentry : deserializer_sentry {
    result_reader_sentry(deserializer& s, std::string const& key, result_kind kind)
        : deserializer_sentry(s, key) {
        uint64_t version; uint32_t stored_kind;
        alps::serialization::deserialize(s, "@version", version);
        alps::serialization::deserialize(s, "@kind", stored_kind);
        if (version != 1 || stored_kind != static_cast<uint32_t>(kind))
            throw std::runtime_error("unsupported ALEA version or statistical kind");
    }
};
}
