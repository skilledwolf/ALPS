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
#include <stdexcept>
#include <boost/filesystem.hpp>
#include <alps/hdf5/archive.hpp>
#include <alps/parameter.h>

using namespace std;

namespace {
void require(bool condition) {
    if (!condition) throw std::runtime_error("legacy Parameters HDF5 round trip failed");
}
}

int main ()
{
    
    if (boost::filesystem::exists("parms.h5") && boost::filesystem::is_regular_file("parms.h5"))
        boost::filesystem::remove("parms.h5");
    
    alps::Parameters p, p2;
    p["a"] = 10;
    p["b"] = "test";
    p["c"] = 10.;
    p["d"] = 5.;
    p["negative"] = -7;
    p["precise"] = "0.12345678901234567";
    
    {
        alps::hdf5::archive ar("parms.h5", "a");
        ar << alps::make_pvp("/parameters", p);
    }
    {
        alps::hdf5::archive ar("parms.h5", "r");
        alps::Parameters pin;
        ar >> alps::make_pvp("/parameters", pin);
        require(ar.is_datatype<int>("/parameters/a"));
        require(ar.is_datatype<std::string>("/parameters/b"));
        require(ar.is_datatype<double>("/parameters/precise"));
        double stored;
        ar.read("/parameters/precise", stored);
        require(pin["a"].get<int>() == 10 && pin["b"] == "test");
        require(pin["c"].get<int>() == 10 && pin["d"].get<int>() == 5);
        require(pin["negative"].get<int>() == -7);
        require(pin["precise"].get<double>() == stored);
        cout << "Reading 1:" << endl << pin;
    }
    
    // "a" is modified from int to double
    // "c" is modified from double to double (but with decimals)
    p2["a"] = 10.5;
    p2["c"] = 5.2;
    {
        alps::hdf5::archive ar("parms.h5", "a");
        ar << alps::make_pvp("/parameters", p2);
    }
    {
        alps::hdf5::archive ar("parms.h5", "r");
        alps::Parameters pin;
        ar >> alps::make_pvp("/parameters", pin);
        require(ar.is_datatype<double>("/parameters/a"));
        require(ar.is_datatype<double>("/parameters/c"));
        require(pin["a"].get<double>() == 10.5 && pin["c"].get<double>() == 5.2);
        require(pin["b"] == "test" && pin["d"].get<int>() == 5);
        cout << "Reading 2:" << endl << pin;
    }
    
    // "d" is modified from double to string
    p2["d"] = "newtype";
    {
        alps::hdf5::archive ar("parms.h5", "a");
        ar << alps::make_pvp("/parameters", p2);
    }
    {
        alps::hdf5::archive ar("parms.h5", "r");
        alps::Parameters pin;
        ar >> alps::make_pvp("/parameters", pin);
        require(ar.is_datatype<std::string>("/parameters/d"));
        require(pin["d"] == "newtype" && pin["b"] == "test");
        require(pin["a"].get<double>() == 10.5 && pin["c"].get<double>() == 5.2);
        cout << "Reading 3:" << endl << pin;
    }
    
    return 0;
}
