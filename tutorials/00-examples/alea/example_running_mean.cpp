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
    const auto values = samples("E");
    std::vector<double> forward(values.size()), reverse(values.size());
    std::partial_sum(values.begin(), values.end(), forward.begin());
    std::partial_sum(values.rbegin(), values.rend(), reverse.rbegin());
    for (size_t i=0; i<values.size(); ++i) {
        forward[i] /= i+1;
        reverse[i] /= values.size()-i;
    }
    std::cout << "The running mean of E is:\n";
    for (double x : forward) std::cout << x << ' ';
    std::cout << "\nThe reverse running mean of E is:\n";
    for (double x : reverse) std::cout << x << ' ';
    std::cout << '\n';
}
