// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#pragma once
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>

namespace native_mc {
// Screen nearly cancelling raw moments in every jackknife evaluation before
// thermal prefactors. This cannot recover precision lost inside signed sums.
// NaN causes the optional estimate to be omitted, preserving raw evidence.
inline auto moment_difference(uint64_t samples, std::size_t bins) {
    double rounding = (double(samples) + bins + 16) * std::numeric_limits<double>::epsilon();
    rounding = rounding < 1 ? 8 * rounding / (1 - rounding) : std::numeric_limits<double>::infinity();
    return [rounding](double a, double b) {
        double scale = std::abs(a) + std::abs(b), delta = a - b;
        return scale && std::abs(delta) <= rounding * scale ? std::numeric_limits<double>::quiet_NaN() : delta;
    };
}
}
