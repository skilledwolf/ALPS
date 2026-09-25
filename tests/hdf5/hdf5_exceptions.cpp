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

#include <alps/hdf5/archive.hpp>

#include <boost/random.hpp>
#include <boost/filesystem.hpp>

#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <stdexcept>

template <typename Exception, typename Operation>
void expect_failure(Operation operation) {
    try {
        operation();
    } catch (Exception const &) {
        return;
    }
    throw std::runtime_error("expected HDF5 operation to fail");
}

int main() {
    std::string const filename = "test_hdf5_exceptions.h5";
    if (boost::filesystem::exists(boost::filesystem::path(filename)))
        boost::filesystem::remove(boost::filesystem::path(filename));
    {
        alps::hdf5::archive oar(filename, "a");
    }
    {
        using namespace alps;
        alps::hdf5::archive iar(filename, "r");
        double test;
        try {
            iar >> make_pvp("/not/existing/path", test);
        } catch (std::exception& ex) {
            std::string str = ex.what();
            std::size_t start = str.find_first_of("\n");
            std::cout << str.substr(0, start) << std::endl;
        }
    }
    {
        alps::hdf5::archive ar(filename, "a");
        ar.create_group("/group");
        ar.write("/dataset", 42);
        for (std::string const parent : {"/group", "/dataset"}) {
            std::string path = parent + "/@values";
            int values[] = {3, 5, 8};
            ar.write(path, values, {3});
            // A rejected partial attribute operation must release every handle.
            expect_failure<std::logic_error>([&] { ar.read(path, values, {2}); });
            expect_failure<std::logic_error>([&] { ar.write(path, values, {3}, {2}); });
            expect_failure<alps::hdf5::archive_error>([&] {
                ar.write(path, static_cast<int const *>(nullptr), {3});
            });
            std::string strings[] = {"one", "two", "three"};
            ar.write(path, strings, {3});
            expect_failure<std::logic_error>([&] { ar.read(path, strings, {2}); });
            // Reuse the same archive after failure, including type/shape changes.
            ar.write(parent + "/@scalar", 17);
            int scalar = 0;
            ar.read(parent + "/@scalar", scalar);
            if (scalar != 17)
                throw std::runtime_error("scalar attribute replacement failed");
            ar.write(path, values, {3});
            int restored[3] = {};
            ar.read(path, restored, {3});
            if (!std::equal(values, values + 3, restored))
                throw std::runtime_error("attribute recovery failed");
        }
        expect_failure<alps::hdf5::path_not_found>([&] { ar.write("/missing/@x", 1); });
        expect_failure<alps::hdf5::archive_error>([&] {
            ar.write("/dataset", static_cast<int const *>(nullptr), {3});
        });
        ar.write("/dataset", 23);
        ar.write("/dataset", std::string("recovered"));
        ar.write("/group", 29);
    } // Closing after all failures must succeed without aborting.
    {
        alps::hdf5::archive ar(filename, "a");
        for (std::string const path : {"/shape", "/dataset/@shape"}) {
            ar.write(path, 17);
            int values[] = {3, 5, 8};
            ar.write(path, values, {3});
            int restored[3] = {};
            ar.read(path, restored, {3});
            if (!std::equal(values, values + 3, restored))
                throw std::runtime_error("scalar to vector replacement failed");
            ar.write(path, 19);
            ar.write(path, values, {0});
            if (!ar.is_null(path))
                throw std::runtime_error("scalar to null replacement failed");
            ar.write(path, std::string("new type"));
        }
    }
    {
        alps::hdf5::archive ar(filename, "r");
        std::string recovered;
        ar.read("/dataset", recovered);
        if (recovered != "recovered")
            throw std::runtime_error("dataset recovery failed");
    }
    boost::filesystem::remove(boost::filesystem::path(filename));
    return 0;
}
