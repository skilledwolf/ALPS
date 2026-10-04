// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
// Compile against the isolated ALPSCore reference SDK, never ALPS headers.
#include <alps/alea.hpp>
#include <alps/alea/hdf5.hpp>
#include <vector>

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
    alps::hdf5::archive archive(argv[1], "w");
    hdf5_serializer bridge(archive, "/results");
    serialize(bridge, "mean", mean.result());
    serialize(bridge, "variance", variance.result());
    serialize(bridge, "covariance", covariance.result());
    serialize(bridge, "autocorr", autocorr.result());
    serialize(bridge, "batch", batch.result());
    serialize(bridge, "circular", circular.result());
    serialize(bridge, "elliptic", elliptic.result());
}
