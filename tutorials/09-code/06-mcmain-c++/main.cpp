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

#include <alps/ngs/api.hpp>
#include <alps/parseargs.hpp>
#include <alps/stop_callback.hpp>
#include "spin_config.hpp"

#include <boost/chrono.hpp>
#include <boost/filesystem/path.hpp>

#include <string>
#include <iostream>
#include <stdexcept>

int main(int argc, char *argv[]) {

    try {
        alps::parseargs options(argc, argv);
        std::string checkpoint_file = options.input_file.substr(0, options.input_file.find_last_of('.')) +  ".clone0.h5";

        alps::params parameters;
        parameters = load_spin_parameters(options.input_file);

        ising_sim sim(parameters);

        if (options.resume)
            sim.load(checkpoint_file);

        if (sim.fraction_completed() < 1.)
            sim.run(alps::stop_callback(options.timelimit));

        sim.save(checkpoint_file);

        auto results = sim.collect_results();

        for (auto const& entry : results)
            std::cout << entry.first << ": " << entry.second << '\n';
        alps::save_results(results, sim.get_parameters(), options.output_file, "/simulation/results");

    } catch (std::exception const & e) {
        std::cerr << "Caught exception: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
