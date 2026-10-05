/*
 * Copyright (C) 1998-2018 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 * For use in publications, see ACKNOWLEDGE.TXT
 */
#include <alps/alea/fisher_f.hpp>
#include <boost/math/distributions/fisher_f.hpp>
#include <cmath>
#include <limits>

namespace alps { namespace alea {
namespace {
double tail(double d1, double d2, double f, bool upper)
{
    if (!std::isfinite(d1) || !std::isfinite(d2) || d1 <= 0 || d2 <= 0 || std::isnan(f))
        return std::numeric_limits<double>::quiet_NaN();
    if (f <= 0) return upper ? 1 : 0;
    if (std::isinf(f)) return upper ? 0 : 1;
    // Reciprocal F swaps its degrees of freedom. Keep Boost's df*f product
    // bounded while retaining tiny upper tails instead of overflowing it.
    if (f > 1) return tail(d2, d1, 1/f, !upper);
    using namespace boost::math::policies;
    using policy_type = policy<domain_error<ignore_error>, evaluation_error<ignore_error>>;
    boost::math::fisher_f_distribution<double, policy_type> dist(d1, d2);
    return upper ? boost::math::cdf(boost::math::complement(dist, f)) : boost::math::cdf(dist, f);
}
}
double fisher_f_distribution::cdf(double f) const { return tail(d1_, d2_, f, false); }
double fisher_f_distribution::ccdf(double f) const { return tail(d1_, d2_, f, true); }
}} /* namespace alps::alea */
