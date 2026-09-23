// SPDX-License-Identifier: MIT
#include <alps/hdf5.hpp>
#include <alps/hdf5/complex.hpp>
#include <alps/hdf5/vector.hpp>
#include <alps/hdf5/valarray.hpp>
#include <alps/hdf5/ublas/vector.hpp>
#include <boost/filesystem/operations.hpp>
#include <complex>
#include <stdexcept>

void require(bool condition) {
    if (!condition) throw std::runtime_error("HDF5 sequence/scalar contract changed");
}

template<class Sequence>
void roundtrip(alps::hdf5::archive& ar, std::string const& path, Sequence const& source) {
    ar[path] << source;
    Sequence restored;
    ar[path] >> restored;
    require(restored.size() == source.size());
    for (std::size_t i = 0; i < source.size(); ++i) require(restored[i] == source[i]);
}

template<template<class...> class Sequence>
void check_sequences(alps::hdf5::archive& ar) {
    Sequence<double> numbers(3);
    numbers[0] = 0.5; numbers[1] = -2; numbers[2] = 4;
    roundtrip(ar, "/values", numbers);
    roundtrip(ar, "/attrs/@values", numbers);
    roundtrip(ar, "/empty", Sequence<double>());
    Sequence<std::complex<double>> complex(2);
    complex[0] = {1, 2}; complex[1] = {-3, 4};
    roundtrip(ar, "/complex", complex);
    roundtrip(ar, "/attrs/@complex", complex);
    Sequence<std::string> strings(2);
    strings[0] = "alpha"; strings[1] = "beta";
    roundtrip(ar, "/strings", strings);
    roundtrip(ar, "/empty-strings", Sequence<std::string>());

    Sequence<std::vector<int>> nested(2);
    nested[0] = {1, 2}; nested[1] = {3, 4};
    roundtrip(ar, "/nested", nested);
    require(ar.is_data("/nested"));
    require(ar.extent("/nested") == std::vector<std::size_t>({2, 2}));
    int flat[4] = {};
    ar.read("/nested", flat, {2, 2}, {0, 0});
    for (int i = 0; i < 4; ++i) require(flat[i] == i + 1);
    nested[1].push_back(5);
    require(!alps::hdf5::is_vectorizable(nested));
    roundtrip(ar, "/nested", nested);
    require(ar.is_group("/nested"));
    // A rectangular write replaces the ragged group's representation.
    nested[1].pop_back();
    roundtrip(ar, "/nested", nested);
    require(ar.is_data("/nested"));

    bool rejected = false;
    try { ar["/values"] >> complex; }
    catch (alps::hdf5::archive_error const&) { rejected = true; }
    require(rejected);
    rejected = false;
    try { alps::hdf5::set_extent(numbers, {}); }
    catch (alps::hdf5::archive_error const&) { rejected = true; }
    require(rejected);
}

template<class T>
void check_scalar(alps::hdf5::archive& ar) {
    for (std::string path : {"/scalar", "/attrs/@scalar"}) {
        ar[path] << T(1);
        T same{};
        double number = 0;
        std::string text;
        ar[path] >> same;
        ar[path] >> number;
        ar[path] >> text;
        require(same == T(1) && number == 1 && std::stold(text) == 1);
    }
}

int main() {
    auto file = boost::filesystem::temp_directory_path() / boost::filesystem::unique_path("alps-codec-%%%%-%%%%.h5");
    {
        alps::hdf5::archive ar(file.string(), "w");
        ar.create_group("/attrs");
        check_sequences<std::vector>(ar);
        check_sequences<boost::numeric::ublas::vector>(ar);
        check_sequences<std::valarray>(ar);
        roundtrip(ar, "/bool", std::vector<bool>{true, false, true});
        roundtrip(ar, "/empty-bool", std::vector<bool>{});
        check_scalar<bool>(ar);
        check_scalar<char>(ar);
        check_scalar<signed char>(ar);
        check_scalar<unsigned char>(ar);
        check_scalar<short>(ar);
        check_scalar<unsigned short>(ar);
        check_scalar<int>(ar);
        check_scalar<unsigned>(ar);
        check_scalar<long>(ar);
        check_scalar<unsigned long>(ar);
        check_scalar<long long>(ar);
        check_scalar<unsigned long long>(ar);
        check_scalar<float>(ar);
        check_scalar<double>(ar);
        check_scalar<long double>(ar);
        for (std::string path : {"/text", "/attrs/@text"}) {
            ar[path] << std::string("123");
            // Repeated attribute reads exercise the shared variable-string lifetime.
            for (int i = 0; i < 100; ++i) {
                int number;
                std::string text;
                ar[path] >> number;
                ar[path] >> text;
                require(number == 123 && text == "123");
            }
            ar[path] << std::string("");
            std::string empty;
            ar[path] >> empty;
            require(empty.empty());
        }
    }
    boost::filesystem::remove(file);
}
