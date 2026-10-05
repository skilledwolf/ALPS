/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2013 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "ising.hpp"


ising_sim::ising_sim(parameters_type const & parms, std::size_t seed_offset)
    : alps::mcbase(parms, seed_offset)
    , length(parameters["L"])
    , sweeps(0)
    , thermalization_sweeps(int(parameters["THERMALIZATION"]))
    , total_sweeps(int(parameters["SWEEPS"]))
    , beta(1. / double(parameters["T"]))
    , spins(length)
{
    for(int i = 0; i < length; ++i)
        spins[i] = (random() < 0.5 ? 1 : -1);
    for (auto const* name : {"Energy", "Magnetization", "Magnetization^2", "Magnetization^4"})
        measurements.emplace(name, std::make_shared<alps::alea::batch_acc<double>>(1, 64));
    measurements.emplace("Correlations", std::make_shared<alps::alea::batch_acc<double>>(length, 64));
}

void ising_sim::update() {
    for (int j = 0; j < length; ++j) {
        using std::exp;
        int i = int(double(length) * random());
        int right = ( i + 1 < length ? i + 1 : 0 );
        int left = ( i - 1 < 0 ? length - 1 : i - 1 );
        double p = exp(-2. * beta * spins[i] * ( spins[right] + spins[left] ));
        if ( p >= 1. || random() < p )
            spins[i] = -spins[i];
    }
}

void ising_sim::measure() {
    sweeps++;
    if (sweeps > thermalization_sweeps) {
        double tmag = 0;
        double ten = 0;
        std::vector<double> corr(length);
        for (int i = 0; i < length; ++i) {
            tmag += spins[i];
            ten += -spins[i] * spins[ i + 1 < length ? i + 1 : 0 ];
            for (int d = 0; d < length; ++d)
                corr[d] += spins[i] * spins[( i + d ) % length ];
        }
        for (double& correlation : corr)
            correlation /= length;
        ten /= length;
        tmag /= length;
        *measurement("Energy") << alps::alea::make_adapter(ten);
        *measurement("Magnetization") << alps::alea::make_adapter(tmag);
        *measurement("Magnetization^2") << alps::alea::make_adapter(tmag * tmag);
        *measurement("Magnetization^4") << alps::alea::make_adapter(tmag * tmag * tmag * tmag);
        *measurement("Correlations") << alps::alea::make_adapter(corr);
    }
}

double ising_sim::fraction_completed() const {
    return (sweeps < thermalization_sweeps ? 0. : ( sweeps - thermalization_sweeps ) / double(total_sweeps));
}

void ising_sim::save(alps::hdf5::archive & ar) const {
    mcbase::save(ar);
    ar["checkpoint/sweeps"] << sweeps;
    ar["checkpoint/spins"] << spins;
}

void ising_sim::load(alps::hdf5::archive & ar) {
    parameters_type restored_parameters;
    ar["/parameters"] >> restored_parameters;
    int restored_length = restored_parameters["L"].as<int>();
    int restored_thermalization = restored_parameters["THERMALIZATION"].as<int>();
    int restored_total = restored_parameters["SWEEPS"].as<int>();
    double restored_beta = 1. / restored_parameters["T"].as<double>();
    int restored_sweeps;
    std::vector<int> restored_spins;
    ar["checkpoint/sweeps"] >> restored_sweeps;
    ar["checkpoint/spins"] >> restored_spins;
    if (restored_length != length || !std::isfinite(restored_beta) || restored_beta <= 0
            || restored_thermalization < 0 || restored_total <= 0 || restored_sweeps < 0
            || std::int64_t(restored_sweeps) > std::int64_t(restored_thermalization) + restored_total
            || restored_spins.size() != std::size_t(restored_length))
        throw std::invalid_argument("invalid Ising checkpoint progress or shape");
    for (auto spin : restored_spins)
        if (spin != -1 && spin != 1)
            throw std::invalid_argument("invalid Ising checkpoint spin");
    mcbase::load(ar);
    thermalization_sweeps = restored_thermalization;
    total_sweeps = restored_total;
    beta = restored_beta;
    sweeps = restored_sweeps;
    spins = std::move(restored_spins);
}
