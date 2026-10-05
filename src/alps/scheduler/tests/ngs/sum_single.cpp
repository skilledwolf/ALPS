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
            *measurement("SValue") << alps::alea::make_adapter(value);
            *measurement("VValue") << alps::alea::make_adapter(std::vector<double>(3, value));
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

    alps::mcoptions options(argc, argv);

    alps::parameters_type<my_sim_type>::type params;
    params = load_sum_parameters(options.input_file);

    my_sim_type my_sim(params); // creat a simulation
    my_sim.run(alps::stop_callback(options.time_limit)); // run the simulation

    alps::results_type<my_sim_type>::type results = collect_results(my_sim); // collect the results

    std::cout << "e^(-x*x): " << std::get<alps::alea::batch_result<double>>(results["SValue"]) << std::endl;
    std::cout << "e^(-x*x): " << std::get<alps::alea::batch_result<double>>(results["VValue"]) << std::endl;
    alps::save_results(results, params, options.output_file, "/simulation/results");
}
