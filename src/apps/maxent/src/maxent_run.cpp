// Copyright (C) 2010 Sebastian Fuchs, Thomas Pruschke, Matthias Troyer.
// Modifications (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#include "maxent.hpp"
#include <alps/hdf5/archive.hpp>
#include <alps/maxent.hpp>
#include <alps/ngs/signal.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
bool alps::solvers::maxent(const params &supplied, const alps::maxent::data &data,
                           const std::string &output_file, int time_limit, bool text_output) {
    if (time_limit < 0)
        throw std::invalid_argument("MaxEnt time limit must be nonnegative");
    const auto parameters = alps::maxent::prepare(supplied, data);
    MaxEntSimulation simulation(parameters, data, output_file, text_output);
    const auto end =
        boost::posix_time::second_clock::local_time() + boost::posix_time::seconds(time_limit);
    alps::ngs::signal signal;
    if (simulation.run([&] {
            return !signal.empty() || boost::posix_time::second_clock::local_time() > end;
        })) {
        alps::hdf5::archive archive(output_file, "a");
        archive["/parameters"] << parameters;
        return true;
    }
    return false;
}
