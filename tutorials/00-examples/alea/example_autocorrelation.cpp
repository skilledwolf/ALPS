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
    const auto values = correlation(samples("m"));
    const auto fit = exponential_fit(values);
    std::cout << "The autocorrelation of m is:\n";
    for (double x : values) std::cout << x << ' ';
    std::cout << "\nThe exponential fit is: " << fit.first << " * exp(" << fit.second << " * t)\n"
              << "The estimated integrated autocorrelation time is: " << integrated_time(values, fit) << '\n';
}
