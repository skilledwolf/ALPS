// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#include <alps/hdf5.hpp>
#include <complex>
#include <filesystem>
#include <stdexcept>
#include <vector>

void require(bool condition) {
    if (!condition)
        throw std::runtime_error("native pointer/pair contract failed");
}

int main() {
    const char *filename = "hdf5_pair.h5";
    {
        alps::hdf5::archive ar(filename, "w");
        const std::vector<std::complex<double>> values{{1, -2}, {3, 4}, {5, -6}, {7, 8}};
        ar << alps::make_pvp("/pointer", values.data(), std::vector<std::size_t>{2, 2});
        require(ar.extent("/pointer") == std::vector<std::size_t>({2, 2}));
        require(ar.is_datatype<std::complex<double>>("/pointer"));
        std::vector<std::complex<double>> actual(values.size());
        ar >> alps::make_pvp("/pointer", actual.data(), std::vector<std::size_t>{2, 2});
        require(actual == values);

        ar << alps::make_pvp("/empty", static_cast<std::complex<double> const *>(nullptr), 0);
        require(ar.extent("/empty") == std::vector<std::size_t>{0} && !ar.is_null("/empty"));
        ar >> alps::make_pvp("/empty", static_cast<std::complex<double> *>(nullptr), 0);
        ar << alps::make_pvp("/empty-matrix", static_cast<double const *>(nullptr),
                            std::vector<std::size_t>{3, 0});
        require(ar.extent("/empty-matrix") == std::vector<std::size_t>({3, 0}));

        std::pair<int, std::complex<double>> pair{42, {3, -4}}, loaded;
        ar["/pair"] << pair;
        ar["/pair"] >> loaded;
        require(loaded == pair && ar.is_scalar("/pair/1"));
        ar["/legacy/first"] << 42;
        ar["/legacy/second"] << 7;
        std::pair<int, int> legacy;
        bool rejected = false;
        try {
            ar["/legacy"] >> legacy;
        } catch (const std::exception &) {
            rejected = true;
        }
        require(rejected);
    }
    std::filesystem::remove(filename);
}
