// Copyright (C) 2026 by the ALPS collaboration
// SPDX-License-Identifier: MIT
// Transitional checks for PR 166 checkpoint fixes. The subsequent native-test
// migration installs the complete ParamsCheckpoint GoogleTest suite.
#include <alps/ngs/params.hpp>
#include <alps/hdf5/vector.hpp>
#include <alps/utility/temporary_filename.hpp>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

void require(bool valid, const char* message) {
    if (!valid) throw std::runtime_error(message);
}
template<class Integer>
void integer(alps::hdf5::archive& archive, Integer value) {
    archive["/integer"] << value;
    archive.set_context("/integer");
    alps::detail::paramvalue scalar;
    scalar.load(archive);
    require(scalar.cast<int>() == static_cast<int>(value), "integer width");
    archive["/integers"] << std::vector<Integer>{0, value};
    archive.set_context("/integers");
    alps::detail::paramvalue array;
    array.load(archive);
    require(array.cast<std::vector<int>>() == std::vector<int>{0, static_cast<int>(value)},
            "integer vector width");
}
int main() {
    const auto filename = alps::temporary_filename("alps-checkpoint-regression-");
    struct Cleanup {
        std::string path;
        ~Cleanup() { std::remove(path.c_str()); }
    } cleanup{filename};
    alps::hdf5::archive archive(filename, "w");
    integer(archive, static_cast<short>(-123));
    integer(archive, static_cast<unsigned short>(123));
    integer(archive, 123U);
    integer(archive, -123L);
    integer(archive, 123UL);
    integer(archive, static_cast<long long>(std::numeric_limits<int>::min()));
    integer(archive, static_cast<unsigned long long>(std::numeric_limits<int>::max()));
    archive["/overflow"] << std::numeric_limits<unsigned long long>::max();
    archive.set_context("/overflow");
    alps::detail::paramvalue existing(23);
    bool overflow = false;
    try { existing.load(archive); }
    catch (std::overflow_error const& error) {
        overflow = std::string(error.what()).find("/overflow") != std::string::npos;
    }
    require(overflow && existing.cast<int>() == 23, "overflow must preserve value and identify path");
    archive["/transaction/a_supported"] << 10;
    archive["/transaction/z_unsupported"] << 1.25F;
    alps::params parameters;
    parameters["keep"] = 42;
    archive.set_context("/transaction");
    bool unsupported = false;
    try { parameters.load(archive); }
    catch (std::runtime_error const& error) {
        unsupported = std::string(error.what()).find("/transaction/z_unsupported") != std::string::npos;
    }
    require(unsupported && parameters.size() == 1 && parameters["keep"].cast<int>() == 42,
            "unsupported load must preserve parameters and identify path");
    archive["/string"] << std::string("dataset");
    archive["/string/@attribute"] << std::string("attribute");
    for (int i = 0; i < 1000; ++i) {
        std::string dataset, attribute;
        archive["/string"] >> dataset;
        archive["/string/@attribute"] >> attribute;
        require(dataset == "dataset" && attribute == "attribute", "scalar string reads");
    }
}
