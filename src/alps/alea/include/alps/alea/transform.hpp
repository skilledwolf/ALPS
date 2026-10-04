/*
 * Copyright (C) 1998-2018 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 * For use in publications, see ACKNOWLEDGE.TXT
 */
#pragma once

#include <alps/alea/core.hpp>
#include <alps/alea/computed.hpp>

#include <alps/alea/mean.hpp>
#include <alps/alea/variance.hpp>
#include <alps/alea/covariance.hpp>
#include <alps/alea/batch.hpp>

#include <alps/alea/propagation.hpp>
#include <alps/alea/convert.hpp>
#include <alps/alea/internal/util.hpp>

#include <random>
#include <complex>
#include <limits>
#include <type_traits>


namespace alps { namespace alea {

namespace internal {
// Two matrix products can round a cancelled variance just below zero. Limit
// correction to their floating-point error bound; invalid covariance survives.
inline bool variance_is_roundoff(double variance, double magnitude, size_t input_size)
{
    return std::isfinite(magnitude) && variance < 0
        && -variance <= 4 * input_size * std::numeric_limits<double>::epsilon() * magnitude;
}
template<class T>
inline void clamp_covariance_roundoff(typename eigen<T>::matrix& covariance,
                              column<double> const& magnitude, size_t input_size)
{
    for (size_t i = 0; i != magnitude.size(); ++i) {
        auto variance = std::real(covariance(i, i));
        if (variance_is_roundoff(variance, magnitude(i), input_size))
            covariance(i, i) = T(0);
    }
}
}

// Treat each complex component as the paired real numerator and denominator.
// Elliptic variance retains their covariance without a dense component matrix.
inline var_result<double> ratio_real_imag(
        var_result<std::complex<double>, elliptic_var> const& in)
{
    internal::check_valid(in);
    if (!in.size() || in.store().data2().size() != in.size()) throw size_mismatch();
    if (!internal::valid_weight_count(in.count(), in.count2()))
        throw std::invalid_argument("Invalid ratio squared-weight count");
    var_result<double> out{var_data<double>(in.size())};
    out.store().count() = in.count();
    out.store().count2() = in.count2();
    out.store().data().fill(std::numeric_limits<double>::quiet_NaN());
    out.store().data2().fill(std::numeric_limits<double>::quiet_NaN());
    if (!in.count()) return out;
    for (size_t i = 0; i != in.size(); ++i) {
        auto numerator = in.mean()(i).real(), denominator = in.mean()(i).imag();
        if (!std::isfinite(numerator) || !std::isfinite(denominator) || denominator == 0)
            throw std::domain_error("Real/imaginary ratio requires finite means and a nonzero denominator");
        auto ratio = numerator / denominator;
        if (!std::isfinite(ratio)) throw std::overflow_error("Real/imaginary ratio is not representable");
        out.store().data()(i) = ratio;
        if (in.observations() <= 1) continue; // No independent error estimate.
        auto const& covariance = in.store().data2()(i);
        if (!isfinite(covariance) || covariance.rere() < 0 || covariance.imim() < 0)
            throw std::domain_error("Invalid real/imaginary ratio covariance");
        auto sign = ratio * (ratio * covariance.imim());
        auto variance = (covariance.rere() - ratio * covariance.reim())
                      + ratio * (ratio * covariance.imim() - covariance.imre());
        auto magnitude = std::abs(covariance.rere())
                       + std::abs(ratio) * (std::abs(covariance.reim()) + std::abs(covariance.imre()))
                       + std::abs(sign);
        if (internal::variance_is_roundoff(variance, magnitude, 2)) variance = 0;
        if (variance < 0) throw std::domain_error("Negative real/imaginary ratio variance");
        out.store().data2()(i) = variance / denominator / denominator;
        if (!std::isfinite(out.store().data2()(i)))
            throw std::overflow_error("Real/imaginary ratio variance is not representable");
    }
    return out;
}

template <typename T, typename InResult>
mean_result<T> transform(no_prop, const transformer<T> &tf, const InResult &in)
{
    static_assert(traits<InResult>::HAVE_MEAN, "result does not have mean");
    static_assert(std::is_same<typename traits<InResult>::value_type, T>::value,
                  "Result and transform types are mismatched");

    if (tf.in_size() != in.size())
        throw size_mismatch();

    mean_result<T> res(mean_data<T>(tf.out_size()));
    res.store().data() = tf(in.mean());
    res.store().count() = in.count();
    return res;
}

template <typename T, typename InResult>
typename std::enable_if<traits<InResult>::HAVE_COV, cov_result<T> >::type transform(linear_prop p, const transformer<T> &tf, const InResult &in)
{
    static_assert(traits<InResult>::HAVE_MEAN, "result does not have mean");
    static_assert(traits<InResult>::HAVE_COV, "result does not have covariance");
    static_assert(std::is_same<typename traits<InResult>::value_type, T>::value,
                  "Result and transform types are mismatched");

    if (tf.in_size() != in.size())
        throw size_mismatch();

    typename eigen<T>::matrix jac = jacobian(tf, in.mean(), p.dx());
    auto covariance = in.cov();

    cov_result<T> res(cov_data<T>(tf.out_size()));
    res.store().data() = tf(in.mean());
    res.store().data2() = jac * covariance * jac.adjoint();
    column<double> magnitude = (jac.cwiseAbs() * covariance.cwiseAbs())
                                     .cwiseProduct(jac.cwiseAbs()).rowwise().sum();
    internal::clamp_covariance_roundoff<T>(res.store().data2(), magnitude, in.size());
    res.store().count() = in.count();
    res.store().count2() = in.count2();
    return res;
}

template <typename T, typename InResult>
typename std::enable_if<!traits<InResult>::HAVE_COV, cov_result<T>>::type transform(linear_prop p, const transformer<T> &tf, const InResult &in)
{
    static_assert(traits<InResult>::HAVE_MEAN, "result does not have mean");
    static_assert(traits<InResult>::HAVE_VAR, "result does not have variance");
    static_assert(std::is_same<typename traits<InResult>::value_type, T>::value,
                  "Result and transform types are mismatched");

    if (tf.in_size() != in.size())
        throw size_mismatch();

    typename eigen<T>::matrix jac = jacobian(tf, in.mean(), p.dx());
    auto variance = in.var();

    cov_result<T> res(cov_data<T>(tf.out_size()));
    res.store().data() = tf(in.mean());
    res.store().data2() = jac * variance.asDiagonal() * jac.adjoint();
    column<double> magnitude = jac.cwiseAbs2() * variance.cwiseAbs();
    internal::clamp_covariance_roundoff<T>(res.store().data2(), magnitude, in.size());
    res.store().count() = in.count();
    res.store().count2() = in.count2();
    return res;
}

template <typename T>
batch_result<T> transform(jackknife_prop, const transformer<T> &tf, const batch_result<T> &in)
{
    if (tf.in_size() != in.size())
        throw size_mismatch();

    batch_result<T> res(jackknife(in.store(), tf));
    return res;
}

}}
