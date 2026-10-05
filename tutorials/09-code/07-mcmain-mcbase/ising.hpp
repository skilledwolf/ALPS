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

#ifndef ALPS_TUTORIAL_ISING_HPP
#define ALPS_TUTORIAL_ISING_HPP

#include <alps/mcbase.hpp>

#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/vector.hpp>

#include <boost/function.hpp>
#include <boost/filesystem/path.hpp>

#include <vector>
#include <string>

class ising_sim : public alps::mcbase {

    public:
    using results_type = std::map<std::string,alps::alea::batch_result<double>>;
    results_type collect_results(result_names_type const& names={}) const {
        return collect_results_as<alps::alea::batch_result<double>>(names);
    }

        
        ising_sim(parameters_type const & parms, std::size_t seed_offset = 0);

        void update() override;
        void measure() override;
        double fraction_completed() const override;

        using alps::mcbase::save;
        void save(alps::hdf5::archive & ar) const override;

        using alps::mcbase::load;
        void load(alps::hdf5::archive & ar) override;

    private:
        
        int length;
        int sweeps;
        int thermalization_sweeps;
        int total_sweeps;
        double beta;
        std::vector<int> spins;
};

#endif
