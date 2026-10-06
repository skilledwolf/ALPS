// Copyright (C) 2010-2012 Lukas Gamper; 2026 ALPS Collaboration.
// SPDX-License-Identifier: MIT
#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/complex.hpp>
#include <alps/hdf5/vector.hpp>
#include <cstdio>
#include <iostream>
#include <stdexcept>

namespace {
void require(bool condition) {
    if (!condition) throw std::runtime_error("Native HDF5 copy contract failed");
}
template<class F> void rejects(F const& operation) {
    try {operation();}
    catch (alps::hdf5::archive_error const&) {return;}
    throw std::runtime_error("Invalid HDF5 copy was accepted");
}
template<class T> T read(alps::hdf5::archive& ar,std::string const& path) {
    T value;ar[path]>>value;return value;
}
}

int main() {
    const char* source_file="test_hdf5_copy.h5";
    const char* target_file="test_hdf5_copy2.h5";
    try {
        const std::vector<std::vector<int>> values{{2},{4,4,4},{6,6,6,6,6},{}};
        const std::vector<std::string> labels{"α","Energy"};
        const std::complex<double> complex(1.,-2.);
        {
            alps::hdf5::archive ar(source_file,"w");
            ar["/dat/vec"]<<values;
            ar["/dat/vec/@foo"]<<10;
            ar["/dat/cpx"]<<complex;
            ar["/dat/empty"]<<std::vector<double>{};
            ar["/dat/labels"]<<labels;
            ar["/dat/@version"]<<2;
            ar["/int"]<<7;
        }
        {
            alps::hdf5::archive source(source_file),target(target_file,"w");
            source.set_context("/dat");target.set_context("/nested");
            source.copy("",target,"copy");
            source.copy("../int",target,"int");
            require(source.get_context()=="/dat" && target.get_context()=="/nested");
            require(read<std::vector<std::vector<int>>>(target,"copy/vec")==values);
            require(read<int>(target,"copy/vec/@foo")==10);
            require(read<std::complex<double>>(target,"copy/cpx")==complex);
            require(target.is_complex("copy/cpx"));
            require(read<std::vector<double>>(target,"copy/empty").empty());
            require(read<std::vector<std::string>>(target,"copy/labels")==labels);
            require(read<int>(target,"copy/@version")==2);
            require(read<int>(target,"int")==7);
            rejects([&] {source.copy("cpx",target,"int");});
            require(read<int>(target,"int")==7);
            rejects([&] {source.copy("missing",target,"missing");});
            rejects([&] {source.copy("@version",target,"attribute");});
            rejects([&] {source.copy("cpx",target,"copy/@attribute");});
            target.close();
            rejects([&] {source.copy("cpx",target,"closed");});
        }
        {
            alps::hdf5::archive source(source_file),target(target_file);
            rejects([&] {source.copy("/dat",target,"/readonly");});
            source.close();
            rejects([&] {source.copy("/dat",target,"/closed");});
        }
        std::remove(source_file);std::remove(target_file);
    } catch (std::exception const& error) {
        std::cerr<<error.what()<<'\n';
        std::remove(source_file);std::remove(target_file);
        return 1;
    }
}
