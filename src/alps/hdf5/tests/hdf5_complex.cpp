// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/array.hpp>
#include <alps/hdf5/complex.hpp>
#include <alps/hdf5/map.hpp>
#include <alps/hdf5/multi_array.hpp>
#include <alps/hdf5/pair.hpp>
#include <alps/hdf5/stdarray.hpp>
#include <alps/hdf5/ublas/matrix.hpp>
#include <alps/hdf5/ublas/vector.hpp>
#include <alps/hdf5/valarray.hpp>
#include <alps/hdf5/vector.hpp>
#include <filesystem>
#include <stdexcept>
#include <type_traits>

void require(bool condition) {
    if (!condition)
        throw std::runtime_error("native complex container contract failed");
}

template <class T> void shape(alps::hdf5::archive &ar, const std::string &path,
                              const T &value, std::vector<std::size_t> expected) {
    ar[path] << value;
    require(ar.extent(path) == expected && !ar.is_null(path));
    require(ar.is_datatype<std::complex<double>>(path));
    require(!ar.is_attribute(path + "/@__complex__"));
}

int main() {
    using complex = std::complex<double>;
    static_assert(std::is_same_v<alps::hdf5::scalar_type<complex>::type, complex>);
    const char *filename = "test_hdf5_complex.h5";
    {
        alps::hdf5::archive ar(filename, "w");
        shape(ar, "/scalar", complex(3, 4), {});
        complex scalar;
        ar["/scalar"] >> scalar;
        require(scalar == complex(3, 4) && ar.is_scalar("/scalar"));
        const std::vector<complex> values{{1, -2}, {3, 4}};
        shape(ar, "/vector", values, {2});
        std::vector<complex> restored;
        ar["/vector"] >> restored;
        require(restored == values);
        shape(ar, "/empty-vector", std::vector<complex>{}, {0});
        ar["/empty-vector"] >> restored;
        require(restored.empty());

        shape(ar, "/array", std::array<complex, 2>{{{1, -2}, {3, 4}}}, {2});
        std::array<complex, 2> array;
        ar["/array"] >> array;
        require(array[0] == values[0] && array[1] == values[1]);
        shape(ar, "/empty-array", std::array<complex, 0>{}, {0});
        std::array<complex, 0> empty_array;
        ar["/empty-array"] >> empty_array;
        shape(ar, "/empty-boost-array", boost::array<complex, 0>{}, {0});
        boost::array<complex, 0> empty_boost_array;
        ar["/empty-boost-array"] >> empty_boost_array;

        std::valarray<complex> valarray(values.data(), values.size()), loaded_valarray;
        shape(ar, "/valarray", valarray, {2});
        ar["/valarray"] >> loaded_valarray;
        require(loaded_valarray.size() == 2 && loaded_valarray[0] == values[0] &&
                loaded_valarray[1] == values[1]);
        shape(ar, "/empty-valarray", std::valarray<complex>{}, {0});
        ar["/empty-valarray"] >> loaded_valarray;
        require(loaded_valarray.size() == 0);

        boost::multi_array<complex, 2> multi(boost::extents[2][3]), loaded_multi;
        for (std::size_t i = 0; i < multi.num_elements(); ++i)
            multi.data()[i] = complex(i, -double(i));
        shape(ar, "/multi", multi, {2, 3});
        ar["/multi"] >> loaded_multi;
        require(loaded_multi.num_elements() == 6 &&
                std::equal(multi.data(), multi.data() + 6, loaded_multi.data()));
        boost::multi_array<complex, 2> empty_multi(boost::extents[0][3]);
        shape(ar, "/empty-multi", empty_multi, {0, 3});
        ar["/empty-multi"] >> loaded_multi;
        require(loaded_multi.shape()[0] == 0 && loaded_multi.shape()[1] == 3);

        boost::numeric::ublas::vector<complex> ublas_vector(0), loaded_ublas_vector(2);
        shape(ar, "/empty-ublas-vector", ublas_vector, {0});
        ar["/empty-ublas-vector"] >> loaded_ublas_vector;
        require(loaded_ublas_vector.size() == 0);
        boost::numeric::ublas::matrix<complex> ublas_matrix(2, 0), loaded_ublas_matrix;
        shape(ar, "/empty-ublas-matrix", ublas_matrix, {2, 0});
        ar["/empty-ublas-matrix"] >> loaded_ublas_matrix;
        require(loaded_ublas_matrix.size1() == 2 && loaded_ublas_matrix.size2() == 0);

        shape(ar, "/empty-fixed-inner", std::vector<std::array<complex, 2>>{}, {0, 2});
        std::vector<std::array<complex, 2>> empty_fixed_inner;
        ar["/empty-fixed-inner"] >> empty_fixed_inner;
        require(empty_fixed_inner.empty());

        ar["/empty-flags"] << std::vector<bool>{};
        require(ar.extent("/empty-flags") == std::vector<std::size_t>{0});
        require(ar.is_datatype<bool>("/empty-flags") && !ar.is_null("/empty-flags"));
        std::vector<bool> flags{true};
        ar["/empty-flags"] >> flags;
        require(flags.empty());
        ar["/empty-custom"] << std::vector<std::pair<int, int>>{};
        require(ar.is_group("/empty-custom") && ar.list_children("/empty-custom").empty());
        ar["/empty-map"] << std::map<std::string, int>{};
        require(ar.is_group("/empty-map") && ar.list_children("/empty-map").empty());
        boost::multi_array<std::vector<complex>, 2> custom(boost::extents[1][2]);
        custom[0][0] = {{1, 2}};
        custom[0][1] = {{3, 4}, {5, 6}};
        ar["/custom-multi"] << custom;
        require(ar.is_group("/custom-multi") && !ar.list_children("/custom-multi").empty());
        custom.resize(boost::extents[0][0]);
        ar["/custom-multi"] << custom;
        require(ar.is_group("/custom-multi") && ar.list_children("/custom-multi").empty());
        ar["/custom-multi"] >> custom;
        require(custom.num_elements() == 0);

        ar["/bad-index/999"] << 4;
        std::vector<complex> guarded{{7, 8}};
        bool rejected = false;
        try {
            ar["/bad-index"] >> guarded;
        } catch (const std::exception &) {
            rejected = true;
        }
        require(rejected && guarded == std::vector<complex>{{7, 8}});

        ar["/kept/child"] << 42;
        const std::string bad_text("a\0b", 3);
        rejected = false;
        try {
            ar["/kept"] << std::vector<std::string>{"valid", bad_text};
        } catch (const std::exception &) {
            rejected = true;
        }
        int child = 0;
        ar["/kept/child"] >> child;
        require(rejected && child == 42);

        ar["/zero-inner"] << std::vector<std::vector<complex>>(2);
        require(ar.extent("/zero-inner") == std::vector<std::size_t>({2, 0}));
        std::vector<std::vector<complex>> zero_inner;
        ar["/zero-inner"] >> zero_inner;
        require(zero_inner.size() == 2 && zero_inner[0].empty() && zero_inner[1].empty());
    }
    std::filesystem::remove(filename);
}
