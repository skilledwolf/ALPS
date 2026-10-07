// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#ifndef ALPS_UTILITY_ROUND_TRIP_HPP
#define ALPS_UTILITY_ROUND_TRIP_HPP

#include <algorithm>
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>
#include <string>
#include <type_traits>

namespace alps {

// The fewest significant digits, at least the stream default of six, that
// read back as x. A value that needs more than six digits is exact at
// digits10 digits unless it needs up to max_digits10.
template <class R> std::string round_trip_string(R x) {
    static_assert(std::is_floating_point_v<R>);
    std::ostringstream out;
    out.imbue(std::locale::classic());
    for (int digits = 6;; digits = std::max(digits + 1, std::numeric_limits<R>::digits10)) {
        out.str({});
        out << std::setprecision(digits) << x;
        std::istringstream in(out.str());
        in.imbue(std::locale::classic());
        R y;
        if (digits >= std::numeric_limits<R>::max_digits10 || ((in >> y) && y == x))
            return out.str();
    }
}

} // namespace alps

#endif
