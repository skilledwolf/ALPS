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

#include <alps/parseargs.hpp>
#include <alps/stop_callback.hpp>
#include <alps/mcmpiadapter.hpp>
#include "spin_config.hpp"

#include <boost/chrono.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem/path.hpp>

#include <string>
#include <iostream>
#include <stdexcept>

int main(int argc, char *argv[]) {

    try {
        boost::mpi::environment env(argc, argv);
        boost::mpi::communicator comm;

        alps::parseargs options(argc, argv);
        std::string checkpoint_file = options.input_file.substr(0, options.input_file.find_last_of('.')) 
                                    +  ".clone" + boost::lexical_cast<std::string>(comm.rank()) + ".h5";

        alps::parameters_type<ising_sim>::type parameters;
        if (comm.rank() > 0)
          /* do nothing*/ ;
        else parameters = load_spin_parameters(options.input_file);

        parameters.broadcast(comm, 0);

        alps::mcmpiadapter<ising_sim> sim(parameters, comm, alps::check_schedule(options.tmin, options.tmax));

        if (options.resume)
            sim.load(checkpoint_file);

        sim.run(alps::stop_callback(comm, options.timelimit));

        sim.save(checkpoint_file);

        using alps::collect_results;
        alps::results_type<ising_sim>::type results = collect_results(sim);

        if (comm.rank() == 0) {
            for (auto const& entry : results)
                std::cout << entry.first << ": " << entry.second << '\n';
            alps::save_results(results, sim.get_parameters(), options.output_file, "/simulation/results");
        }

    } catch (std::exception const & e) {
        std::cerr << "Caught exception: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
