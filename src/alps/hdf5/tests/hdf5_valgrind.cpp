/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2012 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <iostream>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>
#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/vector.hpp>
#include <hdf5.h>

int main() {
    std::string const filename = "test_hdf5_valgrind.h5";
    std::string const overflow(4096, '9');
    std::vector<std::string> const values{"11", overflow, "33", "44", "55", "66"};
    std::vector<std::size_t> const shape{2, 3};
    std::vector<std::size_t> const chunk{1, 2};
    std::vector<double> const doubles(10, 2.);

    try {
        std::remove(filename.c_str());
        ssize_t const objects_before = H5Fget_obj_count(H5F_OBJ_ALL, H5F_OBJ_ALL);
        if (objects_before < 0)
            throw std::runtime_error("could not count HDF5 objects");
        {
            alps::hdf5::archive ar(filename, "w");
            ar.create_group("/group");
            ar.write("/scalar", overflow);
            ar.write("/group/@scalar", overflow);
            ar.write("/scalar/@scalar", overflow);
            ar.write("/vector", values.data(), shape);
            ar.write("/group/@vector", values.data(), shape);
            ar.write("/scalar/@vector", values.data(), shape);
            ar["/doubles"] << doubles;
        }

        auto expect_wrong_type = [&](auto read) {
            auto const before = H5Fget_obj_count(H5F_OBJ_ALL, H5F_OBJ_ALL);
            try {
                read();
            } catch (alps::hdf5::wrong_type const &) {
                if (H5Fget_obj_count(H5F_OBJ_ALL, H5F_OBJ_ALL) != before)
                    throw std::runtime_error("HDF5 object leaked during a rejected read");
                return;
            }
            throw std::runtime_error("string/numeric type mismatch did not throw");
        };

        // Repeat successful and throwing reads so a leak detector can expose
        // unreclaimed HDF5 string buffers, including scalar attributes.
        for (int i = 0; i < 100; ++i) {
            alps::hdf5::archive ar(filename, "r");
            for (char const * path : {"/scalar", "/group/@scalar", "/scalar/@scalar"}) {
                std::string actual;
                ar.read(path, actual);
                if (actual != overflow)
                    throw std::runtime_error("scalar string read changed");
                expect_wrong_type([&] {
                    int value = 0;
                    ar.read(path, value);
                });
            }

            for (char const * path : {"/vector", "/group/@vector", "/scalar/@vector"}) {
                std::vector<std::string> actual(values.size());
                ar.read(path, actual.data(), shape);
                if (actual != values)
                    throw std::runtime_error("vector string read changed");
                expect_wrong_type([&] {
                    std::vector<int> converted(values.size());
                    ar.read(path, converted.data(), shape);
                });
            }

            // The memory buffer has two elements; the file has six. Successful
            // selected reads must reclaim against the memory dataspace.
            std::vector<std::string> selected(2);
            ar.read("/vector", selected.data(), chunk, {1, 1});
            if (selected != std::vector<std::string>{"55", "66"})
                throw std::runtime_error("partial string read changed");
            std::vector<int> converted(2);
            expect_wrong_type([&] {
                ar.read("/vector", converted.data(), chunk, {1, 1});
            });
            expect_wrong_type([&] {
                ar.read("/vector", converted.data(), chunk, {0, 0});
            });

            std::vector<double> actual_doubles;
            ar["/doubles"] >> actual_doubles;
            if (actual_doubles != doubles)
                throw std::runtime_error("numeric vector read changed");
        }

        // Attribute conversions also have to release their parent object IDs
        // on both the group and dataset paths when the conversion throws.
        if (H5Fget_obj_count(H5F_OBJ_ALL, H5F_OBJ_ALL) != objects_before)
            throw std::runtime_error("HDF5 object leaked during a string read");
        if (std::remove(filename.c_str()) != 0)
            throw std::runtime_error("could not remove the test archive");
    } catch (std::exception const & error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
