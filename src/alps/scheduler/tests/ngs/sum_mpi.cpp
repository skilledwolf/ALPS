/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2012 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <alps/ngs.hpp>
#include <alps/mcbase.hpp>
#include <alps/parseargs.hpp>
#include <alps/mcmpiadapter.hpp>
#include <alps/stop_callback.hpp>
#include "sum_config.hpp"

#include <cmath>

// Simulation to measure e^(-x*x)
class my_sim_type : public alps::mcbase {

    public:

        my_sim_type(parameters_type const & params, std::size_t seed_offset = 42)
            : alps::mcbase(params, seed_offset)
            , total_count(params["COUNT"])

        {
            measurements.emplace("SValue", std::make_shared<alps::alea::batch_acc<double>>(1, 64));
            measurements.emplace("VValue", std::make_shared<alps::alea::batch_acc<double>>(3, 64));
        }

        // do the calculation in this function
        void update() override {
            double x = random();
            value = exp(-x * x);
        };

        // do the measurements here
        void measure() override {
            ++count;
            *measurements.at("SValue") << alps::alea::make_adapter(value);
            *measurements.at("VValue") << alps::alea::make_adapter(std::vector<double>(3, value));
        };

        double fraction_completed() const override {
            return count / double(total_count);
        }

    private:
        int count = 0;
        int total_count;
        double value;
};

int main(int argc, char *argv[]) {

    try {

        alps::parseargs options(argc, argv);
        boost::mpi::environment env(argc, argv);
        boost::mpi::communicator c;

        alps::parameters_type<my_sim_type>::type params;
        if (c.rank() > 0)
          /* do nothing*/ ;
        else params = load_sum_parameters(options.input_file);
        broadcast(c, params);

        alps::mcmpiadapter<my_sim_type> my_sim(params, c, alps::check_schedule(options.tmin, options.tmax)); // creat a simulation

        my_sim.run(alps::stop_callback(c, options.timelimit)); // run the simulation

        using alps::collect_results;

        if (c.rank() == 0) { // print the results and save it to hdf5
            alps::results_type<alps::mcmpiadapter<my_sim_type> >::type results = collect_results(my_sim);
            std::cout << "e^(-x*x): " << results["SValue"] << std::endl;
            std::cout << "e^(-x*x): " << results["VValue"] << std::endl;
            alps::save_results(results, params, options.output_file, "/simulation/results");
        } else
            collect_results(my_sim);

    } catch(std::exception & ex) {
        std::cerr << ex.what() << std::endl;
        return -1;
    } catch(...) {
        std::cerr << "Fatal Error: Unknown Exception!\n";
        return -2;
    }
}
