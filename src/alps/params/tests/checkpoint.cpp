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
template <class F> void rejects(F f, const std::string &message = {}) {
    try {
        f();
    } catch (const std::exception &error) {
        require(std::string(error.what()).find(message) != std::string::npos);
        return;
    }
    throw std::runtime_error("invalid checkpoint unexpectedly accepted");
}
int main() {
    const char *file = "params-v2-contract.h5";
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
        p["empty strings"] = std::vector<std::string>{};
        p["empty string"] = std::string{};
        p["strings/with.dots"] = std::vector<std::string>{"a,b", "c", ""};
        p["punctuation/@[]"] = "name is data";
        ar["/parameters"] << p;
        alps::params actual;
        ar["/parameters"] >> actual;
        require(actual == p);
        std::string format;
        ar["/parameters/format"] >> format;
        require(format == "alps.params.v2");
        for (const auto &index : ar.list_children("/parameters/entries")) {
            const auto path = "/parameters/entries/" + index;
            const auto shape = ar.extent(path + "/value");
            require(!ar.is_null(path + "/value"));
            require(shape.size() <= 1);
            require(!ar.is_data(path + "/type"));
            require(!ar.is_attribute(path + "/value/@__complex__"));
            require(!ar.is_attribute(path + "/value/@__alps_type__"));
        }
        p.erase("wide");
        ar["/parameters"] << p;
        ar["/parameters"] >> actual;
        require(actual == p);
        // Unsupported strings must fail before removing an existing checkpoint
        // or writing the earlier entries of a new checkpoint.
        const std::string nul_text("a\0b", 3), nul_key("z\0name", 6);
        auto invalid_text = p, invalid_strings = p, invalid_key = p;
        invalid_text["z"] = nul_text;
        invalid_strings["z"] = std::vector<std::string>{"valid", nul_text};
        invalid_key[nul_key] = "valid";
        for (const auto &invalid : {invalid_text, invalid_strings, invalid_key}) {
            rejects([&] { ar["/parameters"] << invalid; }, "NUL");
            ar["/parameters"] >> actual;
            require(actual == p);
            rejects([&] { ar["/rejected"] << invalid; }, "NUL");
            require(!ar.is_group("/rejected/entries") && !ar.is_data("/rejected/format"));
        }
        alps::params::value_type value("value"), restored_value("value");
        value = "kept";
        ar["/value"] << value;
        value = nul_text;
        rejects([&] { ar["/value"] << value; }, "NUL");
        ar["/value"] >> restored_value;
        require(restored_value.as<std::string>() == "kept");
        rejects([&] { ar["/rejected_value"] << value; }, "NUL");
        require(!ar.is_data("/rejected_value"));
        alps::params empty;
        ar["/empty"] << empty;
        ar["/empty"] >> actual;
        require(actual.empty());
        ar["/legacy/value"] << 42;
        actual["preserve"] = 23;
        rejects([&] { ar["/legacy"] >> actual; });
        require(actual["preserve"].as<int>() == 23);
        ar["/bad"] << p;
        ar["/bad/format"] << std::string("alps.params.v1");
        rejects([&] { ar["/bad"] >> actual; }, "offline conversion");
        require(actual.size() == 1);
        ar["/parameters/entries/1/value"] << short(1);
        rejects([&] { ar["/parameters"] >> actual; });
        require(actual.size() == 1);
        alps::params single;
        single["n"] = 7;
        ar["/bad"] << single;
        ar["/bad/entries/0/value"] << 7.25;
        ar["/bad"] >> actual;
        require(actual["n"].as<double>() == 7.25);
        ar["/bad"] << single;
        ar["/bad/entries/0/value"] << std::vector<std::int64_t>{7};
        ar["/bad"] >> actual;
        require(actual["n"].as<std::vector<std::int64_t>>() == std::vector<std::int64_t>{7});
        ar["/bad"] << single;
        ar["/bad/entries/0/value"] << std::vector<std::vector<double>>{{7}};
        rejects([&] { ar["/bad"] >> actual; });
        require(actual.size() == 1);
        p["unset"];
        rejects([&] { ar["/new"] << p; });
        require(!ar.is_group("/new/entries"));
    }
    std::filesystem::remove(file);
}
