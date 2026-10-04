// SPDX-License-Identifier: MIT
#include <alps/hdf5/archive.hpp>
#include <highfive/highfive.hpp>

#include <array>
#include <complex>
#include <cstdio>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
void require(bool condition) {
    if (!condition) throw std::runtime_error("HighFive interoperability contract failed");
}

template<class Read>
void rejects(Read read) {
    auto const before = H5Fget_obj_count(H5F_OBJ_ALL, H5F_OBJ_ALL);
    try {
        read();
    } catch (alps::hdf5::wrong_type const &) {
        require(H5Fget_obj_count(H5F_OBJ_ALL, H5F_OBJ_ALL) == before);
        return;
    }
    throw std::runtime_error("malformed Boolean enum was accepted");
}

template<class T>
void write_alps(alps::hdf5::archive & ar, std::string const & name) {
    using complex = std::complex<T>;
    std::array<complex, 6> const values{{{1, -2}, {3, -4}, {5, -6}, {7, -8}, {9, -10}, {11, -12}}};
    ar.write(name + "/scalar", values[0]);
    ar.write(name + "/array", values.data(), {2, 3});
    ar.write(name + "/array/@scalar", values[1]);
    ar.write(name + "/array/@array", values.data(), {6});
    ar.write(name + "/empty", static_cast<complex const *>(nullptr), {2, 0});
}

template<class T>
void read_highfive(HighFive::File const & file, std::string const & name) {
    using complex = std::complex<T>;
    auto scalar = file.getDataSet(name + "/scalar");
    require(scalar.getDataType() == HighFive::create_datatype<complex>());
    require(scalar.getSpace().getDimensions().empty());
    require(scalar.read<complex>() == complex(1, -2));
    require(scalar.listAttributeNames().empty());
    auto array = file.getDataSet(name + "/array");
    auto values = array.read<std::vector<std::vector<complex>>>();
    require(values.size() == 2 && values[0].size() == 3 && values[1].size() == 3);
    require(values[0][1] == complex(3, -4) && values[1][2] == complex(11, -12));
    require(array.getAttribute("scalar").read<complex>() == complex(3, -4));
    auto attribute = array.getAttribute("array").read<std::vector<complex>>();
    require(attribute.size() == 6 && attribute.back() == complex(11, -12));
    require(array.listAttributeNames() == std::vector<std::string>({"array", "scalar"}));
    auto empty = file.getDataSet(name + "/empty");
    require(empty.getSpace().getDimensions() == std::vector<std::size_t>({2, 0}));
    require(empty.read<std::vector<std::vector<complex>>>() == std::vector<std::vector<complex>>(2));
}

template<class T>
void write_highfive(HighFive::File & file, std::string const & name) {
    using complex = std::complex<T>;
    file.createDataSet(name + "/scalar", complex(13, -14));
    auto array = file.createDataSet(name + "/array",
        std::vector<std::vector<complex>>{{{1, 2}, {3, 4}, {5, 6}}, {{7, 8}, {9, 10}, {11, 12}}});
    array.createAttribute("scalar", complex(15, -16));
    array.createAttribute("array", std::vector<complex>{{17, 18}, {19, 20}});
    file.createDataSet<complex>(name + "/empty", HighFive::DataSpace({0, 2}));
}

template<class T>
void read_alps(alps::hdf5::archive & ar, std::string const & name) {
    using complex = std::complex<T>;
    complex scalar;
    ar.read(name + "/scalar", scalar);
    require(scalar == complex(13, -14));
    require(ar.is_scalar(name + "/scalar") && ar.extent(name + "/scalar").empty());
    require(ar.is_complex(name + "/scalar") && ar.is_datatype<complex>(name + "/scalar"));
    require(!ar.is_datatype<T>(name + "/scalar"));
    std::array<complex, 6> values{};
    ar.read(name + "/array", values.data(), {2, 3});
    require(values.front() == complex(1, 2) && values.back() == complex(11, 12));
    ar.read(name + "/array/@scalar", scalar);
    require(scalar == complex(15, -16));
    std::array<complex, 2> attribute{};
    ar.read(name + "/array/@array", attribute.data(), {2});
    require(attribute.front() == complex(17, 18) && attribute.back() == complex(19, 20));
    require(!ar.is_null(name + "/empty") && ar.extent(name + "/empty") == std::vector<std::size_t>({0, 2}));
    ar.read(name + "/empty", static_cast<complex *>(nullptr), {0, 2});
}
}

int main() {
    std::string const written = "test_hdf5_highfive_alps.h5", external = "test_hdf5_highfive_external.h5";
    std::remove(written.c_str());
    std::remove(external.c_str());
    try {
        auto const before = H5Fget_obj_count(H5F_OBJ_ALL, H5F_OBJ_ALL);
        {
            alps::hdf5::archive ar(written, "w");
            write_alps<float>(ar, "/float");
            write_alps<double>(ar, "/double");
            write_alps<long double>(ar, "/long_double");
            std::array<bool, 3> const flags{true, false, true};
            ar.write("/flag", true);
            ar.write("/flags", flags.data(), {3});
            ar.write("/flags/@flag", false);
            ar.write("/flags/@flags", flags.data(), {3});
        }
        {
            HighFive::File file(written, HighFive::File::ReadOnly);
            read_highfive<float>(file, "/float");
            read_highfive<double>(file, "/double");
            read_highfive<long double>(file, "/long_double");
            require(file.getDataSet("/flag").read<bool>());
            auto flags = file.getDataSet("/flags");
            require(flags.getDataType().getClass() == HighFive::DataTypeClass::Enum);
            require(flags.read<std::vector<bool>>() == std::vector<bool>({true, false, true}));
            require(!flags.getAttribute("flag").read<bool>());
            require(flags.getAttribute("flags").read<std::vector<bool>>() == std::vector<bool>({true, false, true}));
        }
        {
            HighFive::File file(external, HighFive::File::Overwrite);
            write_highfive<float>(file, "/float");
            write_highfive<double>(file, "/double");
            write_highfive<long double>(file, "/long_double");
            auto flags = file.createDataSet("/flags", std::vector<bool>{false, true, false});
            flags.createAttribute("flag", true);
            // HDF5 permits raw enum bytes which name no enum member. They
            // must never become invalid C++ bool object representations.
            auto type = HighFive::create_datatype<bool>();
            std::array<unsigned char, 3> const invalid{0, 2, 1};
            auto malformed = file.createDataSet("/malformed", HighFive::DataSpace({3}), type);
            malformed.write_raw(invalid.data(), type);
            auto attribute = flags.createAttribute("malformed", HighFive::DataSpace({3}), type);
            attribute.write_raw(invalid.data(), type);
            auto scalar = file.createDataSet("/malformed_scalar", HighFive::DataSpace::Scalar(), type);
            scalar.write_raw(&invalid[1], type);
        }
        {
            alps::hdf5::archive ar(external, "r");
            read_alps<float>(ar, "/float");
            read_alps<double>(ar, "/double");
            read_alps<long double>(ar, "/long_double");
            std::array<bool, 3> flags{};
            ar.read("/flags", flags.data(), {3});
            require(flags == std::array<bool, 3>{false, true, false});
            bool flag = false;
            ar.read("/flags/@flag", flag);
            require(flag);
            for (auto path : {"/malformed", "/flags/@malformed"}) {
                std::array<bool, 3> preserved{true, true, false};
                rejects([&] { ar.read(path, preserved.data(), {3}); });
                require(preserved == std::array<bool, 3>{true, true, false});
            }
            flag = false;
            rejects([&] { ar.read("/malformed_scalar", flag); });
            require(!flag);
        }
        require(H5Fget_obj_count(H5F_OBJ_ALL, H5F_OBJ_ALL) == before);
        require(std::remove(written.c_str()) == 0 && std::remove(external.c_str()) == 0);
    } catch (std::exception const & error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
