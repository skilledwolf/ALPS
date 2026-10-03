// SPDX-License-Identifier: MIT
#include <alps/params.hpp>
#include <alps/hdf5/archive.hpp>
#include <boost/serialization/complex.hpp>
#include <filesystem>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <sstream>
#include <stdexcept>
#include <iostream>

void require(bool condition) {
    if (!condition) throw std::runtime_error("Typed params SDK contract failed");
}

int main() {
    try {
    alps::params parameters;
    parameters["count"] = 3;
    parameters["values"] = std::vector<double>{1., 2.};
    parameters["label"] = "sample";
    require(parameters.size() == 3 && !parameters.exists("absent"));
    require(parameters.value_or("absent",7) == 7 && parameters.find("absent") == parameters.end());
    const alps::params copy(parameters);
    require(copy["count"].as<int>() == 3 && copy.begin()->first == "count");
    bool caught = false;
    try { (void)copy["absent"].as<int>(); }
    catch (const std::runtime_error&) { caught = true; }
    require(caught);

    const char* filename = "params-component-contract.h5";
    {
        alps::hdf5::archive archive(filename, "w");
        archive["/parameters"] << parameters;
    }
    {
        alps::hdf5::archive archive(filename, "r");
        alps::params restored; archive["/parameters"] >> restored;
        require(restored["values"].as<std::vector<double>>() == std::vector<double>({1., 2.}));
        require(restored["label"].as<std::string>() == "sample");
    }
    std::filesystem::remove(filename);
    std::stringstream buffer;
    { boost::archive::text_oarchive archive(buffer); archive << copy; }
    alps::params restored;
    { boost::archive::text_iarchive archive(buffer); archive >> restored; }
    require(restored["count"].as<int>() == 3);
    restored.erase("count");
    require(restored.size() == 2 && !restored.exists("count"));
    std::ostringstream output;
    output << restored;
    require(output.str().find("label = sample") != std::string::npos);
    } catch (const std::exception &error) {
        std::cerr << "Typed params SDK contract: " << error.what() << '\n';
        return 1;
    }
}
