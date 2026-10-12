// SPDX-License-Identifier: MIT
// Bespoke validation and owning adapter for the legacy data reader in
// maxent_parms.cpp (Fuchs, Pruschke and Troyer, 2010). The data layout follows
// ALPS reconciliation 90a36c0d4, src/apps/maxent/src/input.cpp.
#include "data.hpp"
#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/vector.hpp>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <stdexcept>

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::invalid_argument(message);
}
void finite(const std::vector<double>& v) {
    for (double x : v) require(std::isfinite(x), "MaxEnt data must be finite");
}
}
void alps::maxent::validate_data(const data& d) {
    const auto n = d.values.size();
    require(n >= 4 && n <= std::size_t(std::numeric_limits<int>::max()),
            "MaxEnt requires between 4 and INT_MAX data values");
    finite(d.values); finite(d.errors); finite(d.covariance); finite(d.tau);
    require(d.tau.empty() || d.tau.size() == n, "tau must match the data count");
    if (d.covariance.empty()) {
        require(d.errors.size() == n, "errors must match the data count");
        for (double x : d.errors) require(x > 0, "errors must be positive");
    } else {
        require(d.errors.empty(), "supply errors or covariance, not both");
        require(n <= std::numeric_limits<std::size_t>::max() / n &&
                d.covariance.size() == n*n, "covariance must be square");
        double scale = 0;
        for (double x : d.covariance) scale = std::max(scale, std::abs(x));
        const double tolerance = 32 * std::numeric_limits<double>::epsilon() * scale;
        for (std::size_t i = 0; i < n; ++i) {
            require(d.covariance[i*n+i] >= 0, "negative covariance diagonal");
            for (std::size_t j = 0; j < i; ++j)
                require(std::abs(d.covariance[i*n+j] - d.covariance[j*n+i]) <= tolerance,
                        "covariance must be symmetric");
        }
    }
}
int maxent_detail::data_count(const alps::maxent::data& d) {
    alps::maxent::validate_data(d);
    return static_cast<int>(d.values.size());
}
alps::maxent::data maxent_detail::legacy_data(const alps::params& p) {
    const int count = p["NDAT"];
    require(count >= 4, "NDAT must be at least 4");
    const auto n = static_cast<std::size_t>(count);
    const bool covariance = p.defined("COVARIANCE_MATRIX");
    const bool hdf = p.defined("DATA") && (p["DATA_IN_HDF5"] | false);
    alps::maxent::data d;
    d.values.resize(n);
    if (!covariance) d.errors.resize(n);
    if (hdf) {
        alps::hdf5::archive ar(p["DATA"].cast<std::string>(), "r");
        ar["/Data"] >> d.values;
        if (covariance) ar["/Covariance"] >> d.covariance;
        else ar["/Error"] >> d.errors;
        require(d.values.size() == n, "NDAT differs from /Data length");
    } else if (p.defined("DATA")) {
        std::ifstream in(p["DATA"].cast<std::string>());
        require(bool(in), "cannot open MaxEnt data file");
        std::vector<bool> seen(n);
        while (in >> std::ws && !in.eof()) {
            int i; double x, error;
            require(bool(in >> i >> x >> error), "malformed data row (index value error)");
            require(i >= 0 && i < count && !seen[i], "invalid or duplicate data index");
            seen[i] = true; d.values[i] = x;
            if (!covariance) d.errors[i] = error;
        }
        require(std::all_of(seen.begin(), seen.end(), [](bool x) { return x; }), "incomplete data file");
    } else {
        for (int i = 0; i < count; ++i) {
            d.values[i] = p["X_" + std::to_string(i)];
            if (!covariance) d.errors[i] = p["SIGMA_" + std::to_string(i)];
        }
    }
    if (covariance && !hdf) {
        require(n <= std::numeric_limits<std::size_t>::max()/n, "covariance size overflow");
        std::ifstream in(p["COVARIANCE_MATRIX"].cast<std::string>());
        require(bool(in), "cannot open covariance file");
        d.covariance.resize(n*n); std::vector<bool> seen(n*n);
        while (in >> std::ws && !in.eof()) {
            int i, j; double x;
            require(bool(in >> i >> j >> x), "malformed covariance row (i j value)");
            require(i >= 0 && j >= 0 && i < count && j < count, "invalid covariance index");
            const auto k = std::size_t(i)*n+j;
            require(!seen[k], "duplicate covariance index"); seen[k] = true;
            d.covariance[k] = x;
        }
        require(std::all_of(seen.begin(), seen.end(), [](bool x) { return x; }), "incomplete covariance file");
    }
    bool explicit_tau = false;
    for (int i = 0; i < count; ++i) explicit_tau |= p.defined("TAU_" + std::to_string(i));
    if (explicit_tau) for (int i = 0; i < count; ++i) d.tau.push_back(p["TAU_" + std::to_string(i)]);
    alps::maxent::validate_data(d);
    return d;
}
