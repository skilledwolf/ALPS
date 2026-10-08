/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2013 by Lukas Gamper <gamperl@gmail.com>,                  *
 *                              Synge Todo <wistaria@comp-phys.org>                *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <alps/stop_callback.hpp>

namespace alps {

    stop_callback::stop_callback(std::size_t timelimit)
        : limit(timelimit)
        , start(boost::chrono::high_resolution_clock::now())
    {}

    bool stop_callback::operator()() {
        return !signals.empty() || (limit.count() > 0 && boost::chrono::high_resolution_clock::now() > start + limit);
    }

}
