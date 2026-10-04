// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
#include <alps/alea.hpp>
#include <stdexcept>

int main() {
    alps::alea::mean_acc<double> mean;
    mean << alps::alea::make_adapter(1.) << alps::alea::make_adapter(3.);
    if (mean.count() != 2 || mean.result().mean()[0] != 2.)
        throw std::runtime_error("independent ALEA consumer");
}
