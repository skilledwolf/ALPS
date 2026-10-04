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
#include "spin_config.hpp"
#include <alps/alea/hdf5.hpp>
#include <alps/ngs/signal.hpp>
#include <chrono>
#include <filesystem>
#include <iostream>

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: ising timelimit parameter-file\n";
        return EXIT_FAILURE;
    }
    try {
        auto parameters = load_spin_parameters(argv[2]);
        auto basename = std::filesystem::path(argv[2]).replace_extension().string();
        auto checkpoint = basename + ".clone0.h5";
        ising_sim simulation(parameters);
        if (std::filesystem::exists(checkpoint)) simulation.load(checkpoint);
        auto limit = std::chrono::seconds(std::stoul(argv[1]));
        auto started = std::chrono::steady_clock::now();
        alps::ngs::signal signals;
        if (simulation.fraction_completed() < 1)
            simulation.run([&] {
                return !signals.empty() || (limit.count() && std::chrono::steady_clock::now() - started > limit);
            });
        simulation.save(checkpoint);
        auto results = simulation.collect_results();
        alps::hdf5::save_checkpoint(basename + ".out.h5", [&](alps::hdf5::archive& archive) {
            archive["/parameters"] << simulation.get_parameters();
            alps::alea::hdf5_serializer bridge(archive, "/simulation/results");
            for (auto const& entry : results) {
                std::cout << entry.first << ": " << entry.second << '\n';
                serialize(bridge, entry.first, entry.second);
            }
        });
    } catch (std::exception const& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
