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

#include <boost/filesystem.hpp>

#include <string>
#include <stdexcept>

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
        double test = 42.;
        bool caught = false;
        try {
            iar >> make_pvp("/not/existing/path", test);
        } catch (alps::hdf5::path_not_found const& ex) {
            caught = true;
            if (std::string(ex.what()).find("/not/existing/path") == std::string::npos)
                throw std::runtime_error("missing-path diagnostic omits the requested path");
        }
        if (!caught || test != 42.)
            throw std::runtime_error("missing dataset did not preserve the destination and throw path_not_found");
    }
    boost::filesystem::remove(boost::filesystem::path(filename));
    return 0;
}
