// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
#include <alps/alea.hpp>
#include <alps/alea/hdf5.hpp>
#include <stdexcept>
#include <vector>

template<class Accumulator>
void check(alps::alea::hdf5_serializer& reader, char const* name, Accumulator const& expected) {
    decltype(expected.result()) actual;
    deserialize(reader, name, actual);
    if (actual != expected.result()) throw std::runtime_error(name);
}

int main(int argc, char** argv) {
    if (argc != 2) return 1;
    using namespace alps::alea;
    mean_acc<double> mean(3);
    var_acc<double> variance(3, 3);
    cov_acc<double> covariance(3, 3);
    autocorr_acc<double> autocorr(3, 3);
    batch_acc<double> batch(3, 16, 3);
    cov_acc<std::complex<double>, circular_var> circular(2, 3);
    cov_acc<std::complex<double>, elliptic_var> elliptic(2, 3);
    for (int i = 0; i < 67; ++i) {
        std::vector<double> values{double(i), double(i), double(i % 7)};
        auto data = make_adapter(values);
        mean << data; variance << data; covariance << data; autocorr << data; batch << data;
        std::vector<std::complex<double>> z{{double(i), double(i % 7)}, {double(i), -double(i % 7)}};
        circular << make_adapter(z); elliptic << make_adapter(z);
    }
    alps::hdf5::archive archive(argv[1], "r");
    hdf5_serializer reader(archive, "/results");
    check(reader, "mean", mean);
    check(reader, "variance", variance);
    check(reader, "covariance", covariance);
    check(reader, "autocorr", autocorr);
    check(reader, "batch", batch);
    check(reader, "circular", circular);
    check(reader, "elliptic", elliptic);
    archive.close();
}
