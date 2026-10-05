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
    const auto values = samples("m");
    const auto result = statistics(values);
    const auto correlations = correlation(values);
    const double tau = integrated_time(correlations, exponential_fit(correlations));
    const double independent = result.level(0).stderror()(0);
    const double corrected = independent*std::sqrt(1+2*tau);
    std::cout << "The estimated integrated autocorrelation time is: " << tau << '\n'
              << "uncorrelated: " << independent << '\n'
              << "with binning: " << result.stderror()(0) << '\n'
              << "with correlation time: " << corrected << '\n';
    alps::hdf5::archive archive("timeseries.h5", "a");
    archive["/analysis/m/error"] << corrected;
}
