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

#pragma once
#include <alps/alea/autocorr.hpp>
#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/vector.hpp>
#include <unsupported/Eigen/FFT>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <vector>

// These lessons operate on chronological observations, never compressed bins.
inline std::vector<double> samples(std::string const& name) {
    alps::hdf5::archive archive("timeseries.h5");
    std::vector<double> values;
    archive["/samples/"+name] >> values;
    if (values.size()<2 || !std::all_of(values.begin(), values.end(), [](double x) { return std::isfinite(x); }))
        throw std::invalid_argument("At least two finite observations are required");
    std::cout << std::setprecision(17);
    return values;
}
inline auto statistics(std::vector<double> const& values) {
    alps::alea::autocorr_acc<double> accumulator;
    for (double value : values) accumulator << alps::alea::make_adapter(value);
    return accumulator.result();
}
inline std::vector<double> correlation(std::vector<double> const& values) {
    const auto n = values.size();
    const auto result = statistics(values);
    const double mean = result.mean()(0), variance = result.level(0).var()(0);
    if (!(variance>0)) throw std::invalid_argument("Constant observations have undefined correlation");
    std::vector<double> padded(2*n, 0), sums;
    for (size_t i=0; i<n; ++i) padded[i] = values[i]-mean;
    Eigen::FFT<double> fft;
    std::vector<std::complex<double>> spectrum;
    fft.fwd(spectrum, padded);
    for (auto& x : spectrum) x = std::norm(x);
    fft.inv(sums, spectrum);
    std::vector<double> output(n-1);
    for (size_t lag=1; lag<n; ++lag) output[lag-1] = sums[lag]/((n-lag)*variance);
    return output;
}
inline size_t crossing(std::vector<double> const& values, double fraction) {
    return std::find_if(values.begin(), values.end(), [&](double x) { return x<fraction*values.front(); })-values.begin();
}
inline auto exponential_fit(std::vector<double> const& values) {
    const size_t first = crossing(values, .8), stop = crossing(values, .2);
    if (stop<first+2) throw std::invalid_argument("Exponential fit needs at least two positive correlations");
    const double middle = (first+1+stop)/2.;
    double xy=0, xx=0, y=0;
    for (size_t i=first; i<stop; ++i) {
        if (!(values[i]>0)) throw std::invalid_argument("Exponential fit needs positive correlations");
        const double x = i+1-middle, log_y = std::log(values[i]);
        xy += x*log_y; xx += x*x; y += log_y;
    }
    const double slope = xy/xx, amplitude = std::exp(y/(stop-first)-slope*middle);
    if (!(slope<0) || !std::isfinite(amplitude)) throw std::invalid_argument("Tail fit must decay");
    return std::make_pair(amplitude, slope);
}
inline double integrated_time(std::vector<double> const& values, std::pair<double,double> fit) {
    // Include the first crossing, then integrate the exponential tail.
    const size_t count = std::min(values.size(), crossing(values, .2)+1);
    return std::accumulate(values.begin(), values.begin()+count, 0.)
        - fit.first/fit.second*std::exp(fit.second*(count+.5));
}
