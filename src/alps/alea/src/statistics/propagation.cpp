/*
 * Copyright (C) 1998-2018 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 * For use in publications, see ACKNOWLEDGE.TXT
 */
#include <alps/alea/propagation.hpp>

#include <cmath>
#include <limits>
#include <stdexcept>

namespace alps { namespace alea {

template <typename T>
typename eigen<T>::matrix jacobian(const transformer<T> &f, column<T> x, double dx)
{
    size_t in_size = f.in_size();
    size_t out_size = f.out_size();
    if (x.size() != in_size) throw size_mismatch();
    if (!std::isfinite(dx) || dx < 0)
        throw std::invalid_argument("Jacobian step must be finite and nonnegative");
    auto evaluate = [&](column<T> const& point) {
        auto value = f(point);
        if (value.size() != out_size) throw size_mismatch();
        return value;
    };

    typename eigen<T>::matrix result(out_size, in_size);
    for (size_t j = 0; j != in_size; ++j) {
        if (f.is_linear()) {
            column<T> basis = column<T>::Zero(in_size);
            basis(j) = T(1);
            result.col(j) = evaluate(basis);
        } else {
            auto scale = std::abs(x(j));
            auto step = dx ? dx : std::cbrt(std::numeric_limits<double>::epsilon()) * (scale ? scale : 1.);
            column<T> plus = x, minus = x;
            plus(j) += step;
            minus(j) -= step;
            auto width = plus(j) - minus(j);
            if (!std::isfinite(step) || !std::isfinite(std::abs(width)) || width == T(0))
                throw std::invalid_argument("Jacobian step cannot perturb the input");
            result.col(j) = (evaluate(plus) - evaluate(minus)) / width;
        }
    }
    return result;
}

template eigen<double>::matrix jacobian(
            const transformer<double> &, column<double>, double);
template eigen< std::complex<double> >::matrix jacobian(
            const transformer<std::complex<double> > &, column<std::complex<double> >,
            double);


template <typename T>
batch_data<T> jackknife(const batch_data<T> &in, const transformer<T> &tf)
{
    // compute batch sums
    if (tf.in_size() != in.size() || in.count().size() != in.num_batches())
        throw size_mismatch();

    batch_data<T> res(tf.out_size(), in.num_batches());
    column<T> sum_batch = in.batch().rowwise().sum();
    uint64_t sum_count = 0;
    size_t occupied = 0;
    for (size_t i = 0; i != in.num_batches(); ++i) {
        auto count = in.count()(i);
        if (!count && !in.batch().col(i).isZero(0))
            throw std::invalid_argument("Empty jackknife bin has a nonzero sum");
        if (sum_count > std::numeric_limits<uint64_t>::max() - count)
            throw std::overflow_error("Jackknife sample count overflows");
        sum_count += count;
        occupied += count != 0;
    }
    if (occupied < 2)
        throw std::invalid_argument("Jackknife requires at least two occupied bins");
    auto mean_result = tf(column<T>(sum_batch / sum_count));
    if (mean_result.size() != tf.out_size()) throw size_mismatch();

    // compute leave-one-out statistics and transforms
    column<T> leaveout(in.size());
    for (size_t i = 0; i != in.num_batches(); ++i) {
        auto count = in.count()(i);
        if (!count) continue;
        if (tf.is_linear())
            leaveout = in.batch().col(i);
        else
            leaveout = (sum_batch - in.batch().col(i)) / (sum_count - count);
        auto value = tf(leaveout);
        if (value.size() != tf.out_size()) throw size_mismatch();
        // Keep the small bin contribution when the total count is much larger.
        if (tf.is_linear())
            res.batch().col(i) = value;
        else
            res.batch().col(i) = count * mean_result + (sum_count - count) * (mean_result - value);
    }

    res.count() = in.count();

    return res;
}

template batch_data<double> jackknife(const batch_data<double> &in,
                                      const transformer<double> &tf);
template batch_data<std::complex<double> > jackknife(
                                const batch_data<std::complex<double> > &in,
                                const transformer<std::complex<double> > &tf);

}}
