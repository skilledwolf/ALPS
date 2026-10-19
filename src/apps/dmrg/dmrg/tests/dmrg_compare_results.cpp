// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Compares the results of two dmrg runs bit for bit: every dataset below
// /spectrum and /simulation/results (energies, iteration rows, measured
// values and their labels) must have the same shape and the same values.
// Parameters, which differ between an interrupted run and its resumption,
// are not compared.
//
// usage: dmrg_compare_results <first.h5> <second.h5>
#include <alps/hdf5.hpp>

#include <cstddef>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

namespace {

int differences = 0;

void differ(std::string const& path, std::string const& what) {
    std::cerr << "dmrg_compare_results: " << path << ": " << what << '\n';
    ++differences;
}

bool same_bits(std::vector<double> const& a, std::vector<double> const& b) {
    return a.size() == b.size() && (a.empty() || std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0);
}

void compare(alps::hdf5::archive& a, alps::hdf5::archive& b, std::string const& path) {
    if (a.is_group(path)) {
        if (!b.is_group(path)) return differ(path, "missing group");
        auto const children = a.list_children(path);
        if (children != b.list_children(path)) return differ(path, "different members");
        for (auto const& child : children) compare(a, b, path + "/" + child);
        return;
    }
    if (!b.is_data(path)) return differ(path, "missing dataset");
    if (a.extent(path) != b.extent(path)) return differ(path, "different shape");
    if (a.is_datatype<std::string>(path)) {
        std::vector<std::string> x, y;
        a[path] >> x;
        b[path] >> y;
        if (x != y) differ(path, "different strings");
        return;
    }
    if (a.is_scalar(path)) {
        double x = 0., y = 0.;
        a[path] >> x;
        b[path] >> y;
        if (!same_bits({x}, {y})) differ(path, "different value");
        return;
    }
    if (a.dimensions(path) == 1) {
        std::vector<double> x, y;
        a[path] >> x;
        b[path] >> y;
        if (!same_bits(x, y)) differ(path, "different values");
        return;
    }
    if (a.dimensions(path) == 2) {
        std::vector<std::vector<double>> x, y;
        a[path] >> x;
        b[path] >> y;
        if (x.size() != y.size()) return differ(path, "different rows");
        for (std::size_t i = 0; i < x.size(); ++i)
            if (!same_bits(x[i], y[i])) return differ(path, "different values in row " + std::to_string(i));
        return;
    }
    differ(path, "unexpected rank");
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "usage: dmrg_compare_results <first.h5> <second.h5>\n";
        return 2;
    }
    alps::hdf5::archive a(argv[1], "r"), b(argv[2], "r");
    for (std::string const group : {"/spectrum", "/simulation/results"}) {
        if (!a.is_group(group)) differ(group, "missing in the first file");
        else compare(a, b, group);
    }
    if (differences) return 1;
    std::cout << "dmrg_compare_results: identical\n";
    return 0;
}
