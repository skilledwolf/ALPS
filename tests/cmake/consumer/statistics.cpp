// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
#include <alps/alea.hpp>
#include <stdexcept>

int main() {
    alps::alea::mean_acc<double> mean;
    mean << alps::alea::make_adapter(1.) << alps::alea::make_adapter(3.);
    if (mean.count() != 2 || mean.result().mean()[0] != 2.)
        throw std::runtime_error("independent ALEA consumer");
    alps::alea::var_acc<double> variance;
    variance << alps::alea::make_adapter(1.) << alps::alea::make_adapter(3.);
    alps::alea::column<double> expected(1); expected.setZero();
    auto test = alps::alea::test_mean(variance.result(), expected);
    if (std::abs(test.score()-4.) > 1e-12 || test.dist().degrees_of_freedom2() != 1)
        throw std::runtime_error("independent ALEA mean-test consumer");
}
