/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* Copyright (C) 2011-2012 by Lukas Gamper <gamperl@gmail.com>,
*                            Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Maximilian Poprawe <poprawem@ethz.ch>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include "analysis.hpp"

int main() {
    const auto result = statistics(samples("E"));
    const double mean = result.mean()(0);
    std::cout << "The mean of E is: " << mean << '\n';
    alps::hdf5::archive archive("timeseries.h5", "a");
    archive["/analysis/E/mean"] << mean;
}
