// SPDX-License-Identifier: MIT
#include <alps/hdf5.hpp>
#include <alps/hdf5/complex.hpp>
#include <alps/hdf5/vector.hpp>
#include <alps/hdf5/valarray.hpp>
#include <alps/hdf5/ublas/vector.hpp>
#include <alps/hdf5/array.hpp>
#include <alps/hdf5/stdarray.hpp>
#include <alps/hdf5/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
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

template<template<class, std::size_t> class Array>
void check_arrays(alps::hdf5::archive& ar) {
    Array<double, 3> numbers{{0.5, -2, 4}};
    roundtrip(ar, "/array", numbers);
    roundtrip(ar, "/attrs/@array", numbers);
    roundtrip(ar, "/empty-array", Array<int, 0>{});
    Array<std::string, 2> strings{{"alpha", "beta"}};
    roundtrip(ar, "/array-strings", strings);
    Array<Array<std::complex<double>, 2>, 2> complex{};
    complex[0][0] = {1,2}; complex[0][1] = {3,4};
    complex[1][0] = {5,6}; complex[1][1] = {7,8};
    roundtrip(ar, "/array-complex", complex);
    require(ar.extent("/array-complex") == std::vector<std::size_t>({2,2,2}));

    Array<std::vector<int>, 2> nested{{{1,2}, {3,4}}};
    roundtrip(ar, "/array-nested", nested);
    require(ar.is_data("/array-nested"));
    nested[1].push_back(5);
    roundtrip(ar, "/array-nested", nested);
    require(ar.is_group("/array-nested"));
    bool rejected = false;
    try { Array<std::vector<int>, 1> wrong; ar["/array-nested"] >> wrong; }
    catch (alps::hdf5::invalid_path const&) { rejected = true; }
    require(rejected);
    nested[1].pop_back();
    roundtrip(ar, "/array-nested", nested);
    require(ar.is_data("/array-nested"));
    rejected = false;
    try { Array<double, 2> wrong; ar["/array"] >> wrong; }
    catch (alps::hdf5::archive_error const&) { rejected = true; }
    require(rejected);

    // Fixed extent is validated against the current dimension of a hyperslab.
    alps::hdf5::save(ar, "/array-slabs", numbers, {2}, {1}, {0});
    numbers[0] = 9;
    alps::hdf5::save(ar, "/array-slabs", numbers, {2}, {1}, {1});
    Array<double, 3> restored{};
    alps::hdf5::load(ar, "/array-slabs", restored, {1}, {1});
    require(restored == numbers);
    require(alps::hdf5::is_vectorizable(Array<std::vector<int>, 0>{}));
}

void check_tuples(alps::hdf5::archive& ar) {
    auto source = boost::make_tuple(7, std::string("tuple"), std::complex<double>(2, -3), std::vector<int>{4, 5});
    decltype(source) restored;
    ar["/tuple"] << source;
    ar["/tuple"] >> restored;
    require(source == restored && ar.list_children("/tuple").size() == 4);
    require(ar.is_complex("/tuple/2"));
    int first = 0;
    std::string second;
    std::complex<double> third;
    std::vector<int> fourth;
    auto references = boost::tie(first, second, third, fourth);
    ar["/tuple"] >> references;
    require(first == 7 && second == "tuple" && third == std::complex<double>(2, -3));
    require(fourth == std::vector<int>({4, 5}));

    auto ten = boost::make_tuple(0, 1, 2, 3, 4, 5, 6, 7, 8, 9);
    decltype(ten) ten_restored;
    ar["/ten"] << ten;
    ar["/ten"] >> ten_restored;
    require(ten == ten_restored && ar.list_children("/ten").size() == 10);
    auto nested = boost::make_tuple(boost::make_tuple(3, 1.5), boost::make_tuple(-2));
    decltype(nested) nested_restored;
    ar["/tuple-nested"] << nested;
    ar["/tuple-nested"] >> nested_restored;
    require(nested == nested_restored);
    boost::tuple<int, boost::tuples::null_type, int> hole(3, {}, 5), hole_restored;
    ar["/tuple-hole"] << hole;
    ar["/tuple-hole"] >> hole_restored;
    require(boost::get<0>(hole_restored) == 3 && boost::get<2>(hole_restored) == 5);
    require(!ar.is_data("/tuple-hole/1"));
    boost::tuple<> empty;
    ar["/tuple-empty"] << empty;
    ar["/tuple-empty"] >> empty;
    require(!ar.is_group("/tuple-empty") && !ar.is_data("/tuple-empty"));
}

void check_multidimensional(alps::hdf5::archive& ar) {
    boost::multi_array<std::vector<int>, 2> source(boost::extents[2][3]), restored;
    for (int i = 0; i < 6; ++i) source.data()[i] = {i * 10, i * 10 + 1};
    ar["/multi"] << source;
    require(ar.extent("/multi") == std::vector<std::size_t>({2, 3, 2}));
    ar["/multi"] >> restored;
    require(source == restored);
    alps::hdf5::save(ar, "/multi-slabs", source, {2}, {1}, {0});
    source[0][0][0] = 99;
    alps::hdf5::save(ar, "/multi-slabs", source, {2}, {1}, {1});
    alps::hdf5::load(ar, "/multi-slabs", restored, {1}, {1});
    require(source == restored);

    // A raw multidimensional view uses the same coordinates without owning or
    // resizing the outer buffer. Ragged elements retain their group encoding.
    auto view = std::make_pair(source.data(), std::vector<std::size_t>{2, 3});
    auto restored_view = std::make_pair(restored.data(), view.second);
    source[1][0].push_back(72);
    ar["/ragged-view"] << view;
    ar["/ragged-view"] >> restored_view;
    require(source == restored && ar.is_group("/ragged-view"));
    require(ar.is_data("/ragged-view/0/0") && ar.is_data("/ragged-view/1/2"));
    bool rejected = false;
    try { ar["/ragged-multi"] << source; }
    catch (alps::hdf5::wrong_type const&) { rejected = true; }
    require(rejected);

    double values[] = {0, 1, 2, 3, 4, 5}, loaded[6] = {};
    auto numbers = std::make_pair(values, view.second);
    auto loaded_numbers = std::make_pair(loaded, view.second);
    alps::hdf5::save(ar, "/view-slabs", numbers, {2}, {1}, {0});
    values[0] = 9;
    alps::hdf5::save(ar, "/view-slabs", numbers, {2}, {1}, {1});
    alps::hdf5::load(ar, "/view-slabs", loaded_numbers, {1}, {1});
    require(std::equal(std::begin(values), std::end(values), std::begin(loaded)));
    rejected = false;
    try { alps::hdf5::load(ar, "/view-slabs", loaded_numbers, {1, 1, 1, 1}, {}); }
    catch (alps::hdf5::archive_error const&) { rejected = true; }
    require(rejected);

    boost::multi_array<std::vector<std::complex<double>>, 3> complex(boost::extents[2][1][2]), complex_restored;
    for (int i = 0; i < 4; ++i) complex.data()[i] = {{double(i), -2}, {3, double(i)}};
    ar["/multi-complex"] << complex;
    ar["/multi-complex"] >> complex_restored;
    require(complex == complex_restored && ar.is_complex("/multi-complex"));
    require(ar.extent("/multi-complex") == std::vector<std::size_t>({2, 1, 2, 2, 2}));
    alps::multi_array<double, 2> owned(boost::extents[2][3]), owned_restored;
    std::copy(std::begin(values), std::end(values), owned.data());
    ar["/alps-multi"] << owned;
    ar["/alps-multi"] >> owned_restored;
    require(owned == owned_restored);

    boost::multi_array<std::vector<int>, 2> empty(boost::extents[0][3]), empty_restored;
    require(alps::hdf5::is_vectorizable(empty));
    ar["/multi-empty"] << empty;
    ar["/multi-empty"] >> empty_restored;
    require(empty_restored.num_elements() == 0);
    boost::multi_array<std::complex<double>, 2> empty_complex(boost::extents[0][3]), empty_complex_restored;
    ar["/multi-empty-complex"] << empty_complex;
    require(ar.extent("/multi-empty-complex") == std::vector<std::size_t>({0, 3, 2}));
    ar["/multi-empty-complex"] >> empty_complex_restored;
    require(empty_complex_restored.num_elements() == 0);
    rejected = false;
    try { ar["/array"] >> restored; }
    catch (alps::hdf5::archive_error const&) { rejected = true; }
    require(rejected);

    bool boolean = false;
    signed char byte = 0;
    auto bool_view = std::make_pair(&boolean, std::vector<std::size_t>{2, 0});
    auto byte_view = std::make_pair(&byte, bool_view.second);
    ar["/attrs/@empty-bool"] << bool_view;
    ar["/attrs/@empty-byte"] << byte_view;
    require(ar.extent("/attrs/@empty-bool") == bool_view.second);
    require(ar.extent("/attrs/@empty-byte") == byte_view.second);
    auto complex_view = std::make_pair(static_cast<std::complex<double>*>(nullptr), std::vector<std::size_t>{2, 0, 3});
    ar["/attrs/@empty-complex"] << complex_view;
    require(ar.extent("/attrs/@empty-complex") == std::vector<std::size_t>({2, 0, 3, 2}));
    ar["/attrs/@empty-complex"] >> complex_view;
}

int main() {
    auto file = boost::filesystem::temp_directory_path() / boost::filesystem::unique_path("alps-codec-%%%%-%%%%.h5");
    {
        alps::hdf5::archive ar(file.string(), "w");
        ar.create_group("/attrs");
        check_sequences<std::vector>(ar);
        check_sequences<boost::numeric::ublas::vector>(ar);
        check_sequences<std::valarray>(ar);
        check_arrays<std::array>(ar);
        check_arrays<boost::array>(ar);
        check_tuples(ar);
        check_multidimensional(ar);
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
