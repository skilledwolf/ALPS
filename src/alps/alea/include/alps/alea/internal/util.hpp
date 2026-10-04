/*
 * Copyright (C) 1998-2018 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 * For use in publications, see ACKNOWLEDGE.TXT
 */

/**
 * Set of auxiliary processing functions useful for implementations
 */
#pragma once

#include <alps/alea/core.hpp>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace alps { namespace alea { namespace internal {

template <typename Acc>
inline void check_valid(const Acc &acc)
{
    if (!acc.valid())
        throw alps::alea::finalized_accumulator();
}

// Make shape/validity failures collective before exposing any data buffers.
inline reducer_setup check_reduction(const reducer &r, bool valid, size_t size)
{
    auto setup = r.get_setup();
    if (!setup.count || setup.pos >= setup.count)
        throw std::invalid_argument("invalid ALEA reducer setup");
    if (r.get_max(!valid)) throw finalized_accumulator();
    if (r.get_max(!size || size > size_t(std::numeric_limits<int64_t>::max())))
        throw size_mismatch();
    auto dimensions = static_cast<int64_t>(size);
    auto maximum = r.get_max(dimensions);
    auto minimum = -r.get_max(-dimensions);
    if (maximum != minimum) throw size_mismatch();
    return setup;
}

inline bool valid_weight_count(uint64_t count, double count2)
{
    return std::isfinite(count2) && count2 >= 0
        && ((count == 0) == (count2 == 0));
}


template <typename T, typename... Args>
T call_vargs(std::function<T(Args...)> func, const T *args);

template <typename T, typename... Args>
T call_vargs(std::function<T(T, Args...)> func, const T *args)
{
    // use currying to transform multi-argument function to hierarchy
    const T head = *args;
    std::function<T(Args...)> closure =
                [=](Args... tail) -> T { return func(head, tail...); };
    return call_vargs(closure, ++args);
}

template <typename T>
T call_vargs(std::function<T()> func, const T *)
{
    // unravel the currying hierarchy
    return func();
}

template <typename Acc>
typename traits<Acc>::result_type finalize(Acc &acc)
{
    typename traits<Acc>::result_type result;
    acc.finalize_to(result);
    return result;
}

template <typename Acc>
typename traits<Acc>::result_type result(const Acc &acc)
{
    typename traits<Acc>::result_type result;
    Acc copy = acc;
    copy.finalize_to(result);
    return result;
}

}}} /* namespace alps::alea::internal */
