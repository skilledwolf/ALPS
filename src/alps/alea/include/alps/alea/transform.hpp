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

#include <random>
#include <complex>
#include <limits>
#include <type_traits>


namespace alps { namespace alea {

namespace internal {
// Two matrix products can round a cancelled variance just below zero. Limit
// correction to their floating-point error bound; invalid covariance survives.
template<class T>
inline void clamp_covariance_roundoff(typename eigen<T>::matrix& covariance,
                              column<double> const& magnitude, size_t input_size)
{
    auto relative_bound = 4 * input_size * std::numeric_limits<double>::epsilon();
    for (size_t i = 0; i != magnitude.size(); ++i) {
        auto variance = std::real(covariance(i, i));
        if (std::isfinite(magnitude(i)) && variance < 0
                && -variance <= relative_bound * magnitude(i))
            covariance(i, i) = T(0);
    }
}
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
