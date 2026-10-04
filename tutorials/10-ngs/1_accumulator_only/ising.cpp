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

#include <cmath>
#include <sstream>
#include <alps/alea/hdf5.hpp>
#include <alps/hdf5/vector.hpp>

ising_sim::ising_sim(alps::params const & params)
    : parameters(params)
    , random(boost::mt19937((parameters.value_or("SEED", 42))), boost::uniform_real<>())
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
        measurements.emplace(name, alps::alea::batch_acc<double>(1, 64));
    measurements.emplace("Correlations", alps::alea::batch_acc<double>(length, 64));
}

void ising_sim::update() {
    for (int j = 0; j < length; ++j) {
        using std::exp;
        int i = int(double(length) * random());
        int right = ( i + 1 < length ? i + 1 : 0 );
        int left = ( i - 1 < 0 ? length - 1 : i - 1 );
        double p = exp( 2. * beta * spins[i] * ( spins[right] + spins[left] ));
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
        for (auto& value : corr) value /= length;
        ten /= length;
        tmag /= length;
        measurements["Energy"] << alps::alea::make_adapter(ten);
        measurements["Magnetization"] << alps::alea::make_adapter(tmag);
        measurements["Magnetization^2"] << alps::alea::make_adapter(tmag * tmag);
        measurements["Magnetization^4"] << alps::alea::make_adapter(tmag * tmag * tmag * tmag);
        measurements["Correlations"] << alps::alea::make_adapter(corr);
    }
}

double ising_sim::fraction_completed() const {
    return (sweeps < thermalization_sweeps ? 0. : ( sweeps - thermalization_sweeps ) / double(total_sweeps));
}

bool ising_sim::run(std::function<bool()> const& stop_callback) {
    bool stopped = false;
    do {
        update();
        measure();
    } while(!(stopped = stop_callback()) && fraction_completed() < 1.);
    return !stopped;
}

std::map<std::string, alps::alea::batch_result<double>> ising_sim::collect_results() const {
    std::map<std::string, alps::alea::batch_result<double>> results;
    for (auto const& entry : measurements)
        results.emplace(entry.first, entry.second.result());
    return results;
}

void ising_sim::save(std::string const& filename) const {
    alps::hdf5::save_checkpoint(filename, [this](alps::hdf5::archive& archive) {
        save(archive);
    });
}

void ising_sim::load(std::string const& filename) {
    alps::hdf5::archive archive(filename);
    load(archive);
    archive.close();
}

void ising_sim::save(alps::hdf5::archive& archive) const {
    archive["/parameters"] << parameters;
    alps::alea::hdf5_serializer bridge(archive, "/simulation/realizations/0/clones/0/measurements");
    for (auto const& entry : measurements)
        serialize(bridge, entry.first, entry.second);
    archive["/simulation/checkpoint/sweeps"] << sweeps;
    archive["/simulation/checkpoint/spins"] << spins;
    std::ostringstream state;
    state << random.engine();
    archive["/simulation/checkpoint/engine"] << state.str();
}

void ising_sim::load(alps::hdf5::archive& archive) {
    alps::params loaded;
    archive["/parameters"] >> loaded;
    ising_sim restored(loaded);
    alps::alea::hdf5_serializer bridge(archive, "/simulation/realizations/0/clones/0/measurements");
    for (auto& entry : restored.measurements) {
        auto expected = entry.second.size();
        deserialize(bridge, entry.first, entry.second);
        if (entry.second.size() != expected)
            throw std::runtime_error("invalid Ising measurement shape");
    }
    archive["/simulation/checkpoint/sweeps"] >> restored.sweeps;
    archive["/simulation/checkpoint/spins"] >> restored.spins;
    std::string engine;
    archive["/simulation/checkpoint/engine"] >> engine;
    std::istringstream state(engine + " ");
    if (!(state >> restored.random.engine()) || !state.eof()
            || restored.spins.size() != std::size_t(restored.length))
        throw std::runtime_error("invalid Ising checkpoint");
    *this = std::move(restored);
}
