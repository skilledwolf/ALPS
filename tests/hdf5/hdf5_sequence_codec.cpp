// SPDX-License-Identifier: MIT
#include <alps/hdf5.hpp>
#include <alps/hdf5/complex.hpp>
#include <alps/hdf5/vector.hpp>
#include <alps/hdf5/valarray.hpp>
#include <alps/hdf5/ublas/vector.hpp>
#include <alps/hdf5/ublas/matrix.hpp>
#include <alps/hdf5/array.hpp>
#include <alps/hdf5/stdarray.hpp>
#include <alps/hdf5/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/filesystem/operations.hpp>
#include <complex>
#include <limits>
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
        T source[] = {T(0), T(1), T(1), T(0)}, restored[4];
        ar.write(path, source, {2, 2});
        ar.read(path, restored, {2, 2});
        double converted[4];
        ar.read(path, converted, {2, 2});
        for (int i = 0; i < 4; ++i)
            require(restored[i] == source[i] && converted[i] == double(source[i]));
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

template<class Layout>
void check_matrices(alps::hdf5::archive& ar) {
    using matrix = boost::numeric::ublas::matrix<double, Layout>;
    matrix values(2, 3), restored;
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 3; ++j) values(i, j) = 10 * i + j;
    require(alps::hdf5::get_pointer(values) == &values(0, 0));
    require(alps::hdf5::get_pointer(static_cast<matrix const&>(values)) == &values(0, 0));
    for (std::string path : {"/matrix", "/attrs/@matrix"}) {
        ar[path] << values;
        require(ar.extent(path) == std::vector<std::size_t>({2, 3}));
        double raw[6];
        ar.read(path, raw, {2, 3});
        // ALPS stores physical order even for column-major matrices.
        for (int i = 0; i < 6; ++i) require(raw[i] == values.data()[i]);
        ar[path] >> restored;
        require(restored.size1() == 2 && restored.size2() == 3);
        for (std::size_t i = 0; i < 2; ++i)
            for (std::size_t j = 0; j < 3; ++j) require(values(i, j) == restored(i, j));
    }
    alps::hdf5::save(ar, "/matrix-slabs", values, {2}, {1}, {0});
    values(0, 0) = 91;
    alps::hdf5::save(ar, "/matrix-slabs", values, {2}, {1}, {1});
    alps::hdf5::load(ar, "/matrix-slabs", restored, {1}, {1});
    require(restored.size1() == 2 && restored.size2() == 3 && restored(0, 0) == 91);
    for (std::size_t rows : {0, 2}) {
        matrix empty(rows, 0);
        require(alps::hdf5::is_vectorizable(empty) && !alps::hdf5::get_pointer(empty));
        ar["/attrs/@empty-matrix"] << empty;
        ar["/attrs/@empty-matrix"] >> restored;
        require(restored.size1() == rows && restored.size2() == 0);
    }
    using complex_matrix = boost::numeric::ublas::matrix<std::complex<double>, Layout>;
    complex_matrix complex(2, 3), complex_restored;
    for (int i = 0; i < 6; ++i) complex.data()[i] = {double(i), double(-i)};
    ar["/matrix-complex"] << complex;
    ar["/matrix-complex"] >> complex_restored;
    require(ar.extent("/matrix-complex") == std::vector<std::size_t>({2, 3, 2}));
    for (int i = 0; i < 6; ++i) require(complex_restored.data()[i] == complex.data()[i]);
    double guarded[16];
    std::fill(std::begin(guarded), std::end(guarded), 99);
    auto view = std::make_pair(guarded, std::vector<std::size_t>{2, 3});
    auto reject_oversized = [&](auto operation) {
        bool rejected = false;
        try { operation(); }
        catch (alps::hdf5::archive_error const&) { rejected = true; }
        require(rejected);
        for (auto value : guarded) require(value == 99);
    };
    reject_oversized([&] { alps::hdf5::load(ar, "/matrix-complex", view); });
    reject_oversized([&] { alps::hdf5::load(ar, "/matrix-slabs", view, {2}, {0}); });
    reject_oversized([&] { alps::hdf5::save(ar, "/matrix-overflow", view, {2}, {2}, {0}); });
    complex.resize(0, 3);
    ar["/attrs/@empty-complex-matrix"] << complex;
    require(ar.extent("/attrs/@empty-complex-matrix") == std::vector<std::size_t>({0, 3, 2}));
    ar["/attrs/@empty-complex-matrix"] >> complex_restored;
    require(complex_restored.size1() == 0 && complex_restored.size2() == 3);
    bool rejected = false;
    ar["/not-a-matrix"] << std::vector<double>{1, 2, 3};
    try { ar["/not-a-matrix"] >> restored; }
    catch (alps::hdf5::archive_error const&) { rejected = true; }
    require(rejected);
    rejected = false;
    try { ar["/matrix-slabs"] >> restored; }
    catch (alps::hdf5::archive_error const&) { rejected = true; }
    require(rejected); // A rank-three scalar array cannot fit in a matrix.
    boost::numeric::ublas::matrix<std::vector<int>, Layout> unsupported(1, 1);
    rejected = false;
    try { ar["/unsupported-matrix"] << unsupported; }
    catch (alps::hdf5::wrong_type const&) { rejected = true; }
    require(rejected);
}

void check_buffer_conversion(alps::hdf5::archive& ar) {
    std::string empty;
    auto empty_view = std::make_pair(&empty, std::vector<std::size_t>{2, 0});
    ar["/attrs/@empty-string-buffer"] << empty_view;
    require(ar.extent("/attrs/@empty-string-buffer") == empty_view.second);
    ar["/attrs/@empty-string-buffer"] >> empty_view;
    std::string texts[] = {"3", "5", "8", "13", "21", "34"};
    for (std::string path : {"/string-buffer", "/attrs/@string-buffer"}) {
        ar.write(path, texts, {2, 3});
        int numbers[6] = {};
        ar.read(path, numbers, {2, 3});
        require(numbers[0] == 3 && numbers[5] == 34);
    }
    int slab[2] = {};
    ar.read("/string-buffer", slab, {1, 2}, {1, 1});
    require(slab[0] == 21 && slab[1] == 34);
    bool rejected_offset = false;
    std::string invalid_slab[6];
    try { ar.read("/string-buffer", invalid_slab, {2, 3}, {std::numeric_limits<std::size_t>::max(), 0}); }
    catch (alps::hdf5::archive_error const&) { rejected_offset = true; }
    require(rejected_offset);
    std::string invalid[] = {"7", "999999999999999999999999"};
    for (std::string path : {"/invalid-buffer", "/attrs/@invalid-buffer"}) {
        ar.write(path, invalid, {2});
        for (int repeat = 0; repeat < 32; ++repeat) {
            bool rejected = false;
            int numbers[2];
            try { ar.read(path, numbers, {2}); }
            catch (std::out_of_range const&) { rejected = true; }
            require(rejected);
        }
        std::string restored[2];
        ar.read(path, restored, {2});
        require(restored[0] == invalid[0] && restored[1] == invalid[1]);
    }
}

void check_write_replacement(alps::hdf5::archive& ar) {
    for (std::string path : {"/replacement", "/attrs/@replacement"}) {
        ar.write(path, 42);
        roundtrip(ar, path, std::vector<double>{1, 2, 3});
        require(!ar.is_scalar(path) && ar.is_datatype<double>(path));
        roundtrip(ar, path, std::vector<std::string>{"one", "two"});
        require(ar.is_datatype<std::string>(path));
        roundtrip(ar, path, std::vector<std::string>{});
        require(ar.is_null(path));
        roundtrip(ar, path, std::vector<std::string>{"after empty"});
        ar.write(path, std::string("scalar again"));
        std::string text;
        ar.read(path, text);
        require(ar.is_scalar(path) && text == "scalar again");
    }
    // Matching layouts must update in place, preserving dataset attributes.
    int values[] = {1, 2, 3, 4}, patch[] = {8, 9}, restored[4];
    ar.write("/in-place", values, {2, 2});
    ar.write("/in-place/@keep", 7);
    ar.write("/in-place", patch, {2, 2}, {1, 2}, {1, 0});
    ar.read("/in-place", restored, {2, 2});
    require(restored[0] == 1 && restored[1] == 2 && restored[2] == 8 && restored[3] == 9);
    require(ar.is_attribute("/in-place/@keep"));
    ar.write("/in-place", values, {4});
    require(ar.extent("/in-place") == std::vector<std::size_t>{4});
    require(!ar.is_attribute("/in-place/@keep"));

    bool rejected = false;
    try { ar.write("/bad-rank", values, {2, 2}, {2}, {0, 0}); }
    catch (alps::hdf5::archive_error const&) { rejected = true; }
    require(rejected && !ar.is_data("/bad-rank"));
    rejected = false;
    try { ar.write("/attrs/@partial", values, {4}, {2}, {0}); }
    catch (std::logic_error const&) { rejected = true; }
    require(rejected);
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
        check_matrices<boost::numeric::ublas::row_major>(ar);
        check_matrices<boost::numeric::ublas::column_major>(ar);
        check_buffer_conversion(ar);
        check_write_replacement(ar);
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
