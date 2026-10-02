// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/vector.hpp>
#include <alps/params.hpp>
#include <filesystem>
#include <limits>
#include <stdexcept>
void require(bool ok) {
    if (!ok)
        throw std::runtime_error("params checkpoint contract failed");
}
template <class F> void rejects(F f) {
    try {
        f();
    } catch (const std::exception &) {
        return;
    }
    throw std::runtime_error("invalid checkpoint unexpectedly accepted");
}
int main() {
    const char *file = "params-v1-contract.h5";
    {
        alps::hdf5::archive ar(file, "w");
        alps::params p;
        p["wide"] = std::int64_t(9007199254740993LL);
        p["unsigned"] = std::numeric_limits<std::uint64_t>::max();
        p["negative"] = -7;
        p["bool"] = true;
        p["complex"] = std::complex<double>(1, -2);
        p["complex vector"] = std::vector<std::complex<double>>{{1, 2}, {3, -4}};
        p["booleans"] = std::vector<bool>{true, false};
        p["empty booleans"] = std::vector<bool>{};
        p["empty int"] = std::vector<std::int64_t>{};
        p["empty real"] = std::vector<double>{};
        p["empty complex"] = std::vector<std::complex<double>>{};
        p["strings/with.dots"] = std::vector<std::string>{"a,b", "c", ""};
        ar["/parameters"] << p;
        alps::params actual;
        ar["/parameters"] >> actual;
        require(actual == p);
        p.erase("wide");
        ar["/parameters"] << p;
        ar["/parameters"] >> actual;
        require(actual == p);
        alps::params empty;
        ar["/empty"] << empty;
        ar["/empty"] >> actual;
        require(actual.empty());
        ar["/legacy/value"] << 42;
        actual["preserve"] = 23;
        rejects([&] { ar["/legacy"] >> actual; });
        require(actual["preserve"].as<int>() == 23);
        ar["/parameters/entries/1/type"] << std::string("unknown");
        rejects([&] { ar["/parameters"] >> actual; });
        require(actual.size() == 1);
        alps::params single;
        single["n"] = 7;
        ar["/bad"] << single;
        ar["/bad/entries/0/value"] << 7.25;
        rejects([&] { ar["/bad"] >> actual; });
        require(actual.size() == 1);
        ar["/bad"] << single;
        ar["/bad/entries/0/value"] << std::vector<std::int64_t>{7};
        rejects([&] { ar["/bad"] >> actual; });
        require(actual.size() == 1);
        p["unset"];
        rejects([&] { ar["/new"] << p; });
        require(!ar.is_group("/new/entries"));
    }
    std::filesystem::remove(file);
}
