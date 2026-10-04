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

#pragma once

#include <alps/alea.hpp>
#include <alps/alea/checkpoint.hpp>
#include <alps/hdf5/archive.hpp>
#include <alps/params.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <functional>
#include <map>
#include <string>
#include <vector>

class ising_sim {
public:
    explicit ising_sim(alps::params const& parameters);
    void update();
    void measure();
    double fraction_completed() const;
    alps::params const& get_parameters() const { return parameters; }
    bool run(std::function<bool()> const& stop);
    std::map<std::string, alps::alea::batch_result<double>> collect_results() const;
    void save(std::string const& filename) const;
    void load(std::string const& filename);
    void save(alps::hdf5::archive& archive) const;
    void load(alps::hdf5::archive& archive);

private:
    alps::params parameters;
    boost::variate_generator<boost::mt19937, boost::uniform_real<>> random;
    std::map<std::string, alps::alea::batch_acc<double>> measurements;
    int length, sweeps, thermalization_sweeps, total_sweeps;
    double beta;
    std::vector<int> spins;
};
