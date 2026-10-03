// SPDX-License-Identifier: MIT
#include <alps/hdf5/archive.hpp>

#include <hdf5.h>

#include <array>
#include <cstdio>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
void require(bool condition) {
    if (!condition)
        throw std::runtime_error("HDF5 read contract failed");
}

template <typename Exception, typename Read>
void rejects(Read read) {
    try {
        read();
    } catch (Exception const &) {
        return;
    }
    throw std::runtime_error("HDF5 read did not reject invalid input");
}

template <typename T>
void numeric_reads(alps::hdf5::archive & ar) {
    std::array<T, 6> const values{0, 1, 0, 1, 1, 0};
    for (auto path : {"/scalar", "/group/@scalar", "/array/@scalar"}) {
        ar.write(path, T(1));
        double number = -1;
        ar.read(path, number);
        require(number == 1);
        std::string text;
        ar.read(path, text);
        require(std::stod(text) == 1);
    }
    for (auto path : {"/array", "/group/@array", "/scalar/@array"}) {
        ar.write(path, values.data(), {2, 3});
        std::array<double, 6> numbers;
        ar.read(path, numbers.data(), {2, 3});
        require(numbers == std::array<double, 6>{0, 1, 0, 1, 1, 0});
        std::array<std::string, 6> texts;
        ar.read(path, texts.data(), {2, 3});
        for (std::size_t i = 0; i < texts.size(); ++i)
            require(std::stod(texts[i]) == numbers[i]);
    }
    std::array<double, 2> selected;
    ar.read("/array", selected.data(), {2, 1}, {0, 1});
    require(selected == std::array<double, 2>{1, 1});
    ar.read("/array", selected.data(), {1, 2}, {1, 1});
    require(selected == std::array<double, 2>{1, 0});
}

void fixed_strings(std::string const & filename) {
    // External fixed-width strings exercise a representation ALPS does not write.
    hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    hid_t type = H5Tcopy(H5T_C_S1);
    require(file >= 0 && type >= 0 && H5Tset_size(type, 4) >= 0);
    hid_t space = H5Screate(H5S_SCALAR);
    hid_t data = H5Dcreate2(file, "/fixed", type, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hid_t attr = H5Acreate2(data, "fixed", type, space, H5P_DEFAULT, H5P_DEFAULT);
    require(data >= 0 && attr >= 0);
    require(H5Dwrite(data, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, "1234") >= 0);
    require(H5Awrite(attr, type, "1234") >= 0);
    require(H5Aclose(attr) >= 0 && H5Dclose(data) >= 0);
    hsize_t const extent = 2;
    require(H5Sclose(space) >= 0);
    space = H5Screate_simple(1, &extent, nullptr);
    data = H5Dcreate2(file, "/fixed_array", type, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    attr = H5Acreate2(data, "fixed", type, space, H5P_DEFAULT, H5P_DEFAULT);
    require(data >= 0 && attr >= 0);
    require(H5Dwrite(data, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, "12345678") >= 0);
    require(H5Awrite(attr, type, "12345678") >= 0);
    require(H5Aclose(attr) >= 0 && H5Dclose(data) >= 0);
    require(H5Sclose(space) >= 0 && H5Tclose(type) >= 0 && H5Fclose(file) >= 0);

    alps::hdf5::archive ar(filename, "r");
    for (auto path : {"/fixed", "/fixed/@fixed"}) {
        int number = 0;
        ar.read(path, number);
        require(number == 1234);
        std::string text;
        ar.read(path, text);
        require(text == std::string("1234\0", 5));
    }
    for (auto path : {"/fixed_array", "/fixed_array/@fixed"}) {
        std::array<std::string, 2> texts;
        rejects<std::logic_error>([&] { ar.read(path, texts.data(), {2}); });
    }
}
}

int main() {
    std::string const filename = "test_hdf5_read.h5";
    std::remove(filename.c_str());
    try {
        auto const objects_before = H5Fget_obj_count(H5F_OBJ_ALL, H5F_OBJ_ALL);
        {
            alps::hdf5::archive ar(filename, "w");
            ar.create_group("/group");
            int initial[6] = {};
            ar.write("/array", initial, {2, 3});
            numeric_reads<char>(ar);
            numeric_reads<signed char>(ar);
            numeric_reads<unsigned char>(ar);
            numeric_reads<short>(ar);
            numeric_reads<unsigned short>(ar);
            numeric_reads<int>(ar);
            numeric_reads<unsigned>(ar);
            numeric_reads<long>(ar);
            numeric_reads<unsigned long>(ar);
            numeric_reads<long long>(ar);
            numeric_reads<unsigned long long>(ar);
            numeric_reads<float>(ar);
            numeric_reads<double>(ar);
            numeric_reads<long double>(ar);
            numeric_reads<bool>(ar);

            int value = 42;
            std::array<int, 6> values{};
            rejects<alps::hdf5::path_not_found>([&] { ar.read("/missing", value); });
            rejects<alps::hdf5::path_not_found>([&] { ar.read("/group/@missing", value); });
            rejects<alps::hdf5::wrong_type>([&] { ar.read("/array", value); });
            rejects<alps::hdf5::wrong_type>([&] { ar.read("/group/@array", value); });
            rejects<alps::hdf5::archive_error>([&] { ar.read("/scalar", values.data(), {}); });
            rejects<alps::hdf5::archive_error>([&] { ar.read("/group/@scalar", values.data(), {}); });
            rejects<alps::hdf5::archive_error>([&] { ar.read("/array", values.data(), {6}); });
            rejects<alps::hdf5::archive_error>([&] { ar.read("/array", values.data(), {2, 3}, {1, 0}); });
            auto const maximum = std::numeric_limits<std::size_t>::max();
            rejects<alps::hdf5::archive_error>([&] { ar.read("/array", values.data(), {1, 1}, {maximum, 0}); });
            rejects<alps::hdf5::archive_error>([&] { ar.read("/array", values.data(), {maximum, 1}, {1, 0}); });
            require(values == std::array<int, 6>{});
            rejects<alps::hdf5::archive_error>([&] { ar.read("/array", values.data(), {0, 3}); });
            rejects<std::logic_error>([&] { ar.read("/group/@array", values.data(), {1, 3}); });
            ar.write("/empty", static_cast<int const *>(nullptr), {});
            ar.read("/empty", &value, {0}, {0});
            require(value == 42);
            ar.set_context("/group");
            ar.read("@scalar", value);
            require(value == 1);
            ar.close();
            rejects<alps::hdf5::archive_closed>([&] { ar.read("/scalar", value); });
            rejects<alps::hdf5::archive_closed>([&] { ar.read("/array", values.data(), {2, 3}); });
        }
        fixed_strings(filename);
        require(H5Fget_obj_count(H5F_OBJ_ALL, H5F_OBJ_ALL) == objects_before);
        require(std::remove(filename.c_str()) == 0);
    } catch (std::exception const & error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
