// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#include <alps/params.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/complex.hpp>
#include <limits>
#include <sstream>
#include <stdexcept>
using alps::params;
void require(bool ok) {
    if (!ok)
        throw std::runtime_error("params contract failed");
}
template <class F> void rejects(F f) {
    try {
        f();
    } catch (const alps::params_ns::exception::exception_base &) {
        return;
    }
    throw std::runtime_error("parameter conversion unexpectedly accepted");
}
int main() {
    params p;
    const auto &read = p;
    require(!p.exists("missing") && p.value_or("missing", 7) == 7 && p.empty());
    rejects([&] { read["missing"].as<int>(); });
    p["unset"];
    require(!p.exists("unset") && p.size() == 1);
    rejects([&] { p["unset"].as<int>(); });
    p.erase("unset");
    p.erase("missing");
    p["z"] = std::int64_t(9007199254740993LL);
    p["a"] = true;
    require(p.begin()->first == "a" && p["z"].as<std::int64_t>() == 9007199254740993LL);
    rejects([&] { p["z"].as<int>(); });
    rejects([&] { p["z"].as<double>(); });
    rejects([&] { p["a"].as<int>(); });
    p["negative"] = -1;
    rejects([&] { p["negative"].as<unsigned>(); });
    p["unsigned"] = std::numeric_limits<std::uint64_t>::max();
    require(p["unsigned"].as<std::uint64_t>() == std::numeric_limits<std::uint64_t>::max());
    rejects([&] { p["unsigned"].as<std::int64_t>(); });
    p["small"] = 4;
    require(p["small"].as<double>() == 4.);
    p["real"] = 4.5;
    rejects([&] { p["real"].as<int>(); });
    p["text"] = "4";
    rejects([&] { p["text"].as<int>(); });
    p["v"] = std::vector<int>{1, 2};
    require(p["v"].as<std::vector<double>>() == std::vector<double>({1, 2}));
    p["complex"] = std::complex<double>(1, 2);
    require(p["complex"].as<std::complex<double>>() == std::complex<double>(1, 2));
    rejects([&] { p["complex"].as<double>(); });
    auto copy = p;
    copy["small"] = 5;
    require(p["small"].as<int>() == 4);
    std::stringstream buffer;
    {
        boost::archive::text_oarchive ar(buffer);
        ar << p;
    }
    params restored;
    {
        boost::archive::text_iarchive ar(buffer);
        ar >> restored;
    }
    require(restored == p);
}
