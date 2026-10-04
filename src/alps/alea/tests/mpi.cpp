// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#include <alps/alea.hpp>
#include <alps/alea/mpi.hpp>
#include <alps/alea/transform.hpp>
#include <cmath>
#include <complex>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {
namespace aa = alps::alea;
void require(bool condition, char const* message) {
    int local = condition, agreed;
    MPI_Allreduce(&local, &agreed, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if (!agreed) throw std::runtime_error(message);
}
template<class F> void rejects(F operation) {
    bool failed = false;
    try { operation(); } catch (std::exception const&) { failed = true; }
    require(failed, "invalid collective reduction was accepted");
}
template<class R> void rejects_weight_overflow(R result, aa::reducer const& reducer) {
    result.store().count2() = std::numeric_limits<double>::max();
    auto original = result;
    rejects([&] { result.reduce(reducer); });
    require(result == original, "failed collective reduction changed the original result");
}
void empty_autocorrelation(int rank, aa::reducer const& reducer) {
    constexpr size_t samples = 8195;
    auto value = [](size_t i) { return double(int((i/16)%13)-6); };
    aa::autocorr_acc<double> populated;
    for (size_t i=0; i<samples; ++i) populated << value(i);
    auto original = populated.result();
    auto selected = original.find_level(1024);
    require(selected > 0, "empty-rank fixture needs coarse autocorrelation estimates");
    struct oracle { double mean=0., variance=0., count2=0.; };
    auto raw = [&](size_t width) {
        std::vector<double> sums, weights;
        oracle result;
        for (size_t i=0; i<samples; ++i) {
            if (i%width == 0) { sums.push_back(0.); weights.push_back(0.); }
            sums.back() += value(i); ++weights.back();
            result.mean += value(i);
        }
        result.mean /= samples;
        for (size_t i=0; i<sums.size(); ++i) {
            result.count2 += weights[i]*weights[i];
            auto deviation = sums[i]/weights[i]-result.mean;
            result.variance += weights[i]*deviation*deviation;
        }
        result.variance /= samples-result.count2/samples;
        return result;
    };
    auto base = raw(1), coarse = raw(size_t(1)<<selected);
    auto error = std::sqrt(coarse.variance*coarse.count2/(double(samples)*samples));
    auto tau = .5*(coarse.count2/samples)*coarse.variance/base.variance-.5;
    for (int active : {0,1}) {
        aa::autocorr_acc<double> empty(1,9,3);
        auto result = rank == active ? original : empty.result();
        result.reduce(reducer);
        require(result.valid() == (rank == 1), "empty-rank result survived on wrong rank");
        require(rank != 1 || (result.count() == samples && result.nlevel() == original.nlevel()
                 && result.count2() == coarse.count2
                 && std::abs(result.mean()(0)-coarse.mean) < 1e-11
                 && std::abs(result.var()(0)-coarse.variance) < 1e-10
                 && std::abs(result.stderror()(0)-error) < 1e-11
                 && std::abs(result.tau()(0)-tau) < 1e-10),
                "empty rank changed autocorrelation estimates from the raw-bin oracle");
    }
    aa::autocorr_acc<double> empty;
    auto result = empty.result();
    result.reduce(reducer);
    require(result.valid() == (rank == 1), "all-empty result survived on wrong rank");
    require(rank != 1 || (result.count() == 0 && result.nlevel() == 1
                         && result.level(0).count2() == 0.),
            "all-empty autocorrelation reduction invented samples");
    result = empty.result();
    if (!rank) result.level(0).store().count2() = 1.;
    rejects([&] { result.reduce(reducer); });
    require(result.valid() && result.count() == 0 && result.level(0).count2() == (rank ? 0. : 1.),
            "malformed empty-rank rejection changed the original result");
}
void elliptic_signed_ratio(int rank, aa::reducer const& reducer) {
    aa::var_acc<std::complex<double>,aa::elliptic_var> joint(2,2);
    for (int i=0; i<(rank ? 5 : 3); ++i) {
        auto sign = rank && i>=3 ? -1. : 1.;
        auto numerator = sign*(rank ? 5. : 2.);
        joint << std::vector<std::complex<double>>{{numerator,sign},{-2*numerator,sign}};
    }
    auto pooled = joint.result();
    pooled.reduce(reducer);
    require(pooled.valid() == (rank == 1), "elliptic ratio raw state survived on wrong rank");
    bool ok = true;
    if (rank == 1) try {
        auto original = pooled;
        auto result = aa::ratio_real_imag(pooled);
        // Original independent bins: (weight, numerator sum, sign sum).
        double bins[][3]{{2.,4.,2.},{1.,2.,1.},{2.,10.,2.},{2.,0.,0.},{1.,-5.,-1.}};
        double count=0., count2=0., numerator=0., denominator=0.;
        for (auto const& bin : bins) {
            count += bin[0]; count2 += bin[0]*bin[0];
            numerator += bin[1]; denominator += bin[2];
        }
        auto mean = numerator/denominator;
        double variance = 0.;
        for (auto const& bin : bins) {
            auto residual = (bin[1]/bin[0]-mean*bin[2]/bin[0])/(denominator/count);
            variance += bin[0]*residual*residual;
        }
        variance /= count-count2/count;
        auto error = std::sqrt(variance*count2/(count*count));
        ok = pooled == original && result.count() == count && result.count2() == count2
          && result.mean()(0) == mean && result.mean()(1) == -2*mean
          && std::abs(result.var()(0)-variance) < 1e-12
          && std::abs(result.stderror()(0)-error) < 1e-12
          && std::abs(result.stderror()(1)-2*error) < 1e-12;
    } catch (...) { ok = false; }
    require(ok, "elliptic signed ratio lost pooled means, cross-sign covariance or partial-bin weights");
}
void contract(int rank) {
    aa::mpi_reducer reducer(MPI_COMM_WORLD, 1);
    auto setup = reducer.get_setup();
    require(setup.pos == size_t(rank) && setup.count == 2 && setup.have_result == (rank == 1),
            "wrong reduction setup");
    int64_t wide = (int64_t(1) << 33) + rank;
    require(reducer.get_max(wide) == (int64_t(1) << 33) + 1, "get_max narrowed int64_t");
    reducer.reduce(aa::view<int64_t>(&wide, 1));
    require(rank != 1 || wide == (int64_t(1) << 34) + 1, "signed reduction narrowed int64_t");
    uint64_t weight = rank ? 7 : (uint64_t(1) << 63) + 5;
    reducer.reduce(aa::view<uint64_t>(&weight, 1));
    require(rank != 1 || weight == (uint64_t(1) << 63) + 12, "reduction lost full-width uint64_t");
    int32_t small = rank + 1;
    reducer.reduce(aa::view<int32_t>(&small, 1));
    require(rank != 1 || small == 3, "wrong int32_t reduction");
    std::complex<double> complex(rank + 1., 2. * rank + 1.);
    reducer.reduce(aa::view<std::complex<double>>(&complex, 1));
    require(rank != 1 || complex == std::complex<double>(3., 4.), "wrong complex reduction");
    reducer.reduce(aa::view<double>(nullptr, 0));
    rejects([&] { reducer.reduce(aa::view<double>(nullptr, size_t(INT_MAX) + 1)); });
    auto maximum = std::numeric_limits<size_t>::max();
    rejects([&] {
        reducer.reduce(aa::view<std::complex<double>>(nullptr, rank ? maximum / 2 + 1 : 0));
    });
    rejects([&] {
        reducer.reduce(aa::view<aa::complex_op<double>>(nullptr, rank ? maximum / 4 + 1 : 0));
    });
    double value = 1.;
    rejects([&] { reducer.reduce(aa::view<double>(&value, rank)); });
    rejects([&] { aa::mpi_reducer different_roots(MPI_COMM_WORLD, rank); });
    rejects([&] { aa::mpi_reducer invalid_root(MPI_COMM_WORLD, 2); });

    aa::var_acc<double> constant_variance(1);
    aa::cov_acc<double> constant_covariance(1);
    for (int i = 0; i < 2; ++i) {
        constant_variance << 1.;
        constant_covariance << 1.;
    }
    rejects_weight_overflow(constant_variance.result(), reducer);
    rejects_weight_overflow(constant_covariance.result(), reducer);

    aa::mean_acc<double> mean(2);
    aa::var_acc<double> variance(2);
    aa::cov_acc<double> covariance(2);
    aa::batch_acc<double> batch(1, 8, 2);
    aa::autocorr_acc<double> autocorr(1, 3);
    for (int i = rank ? 4 : 1; i <= (rank ? 10 : 3); ++i) {
        std::vector<double> sample{double(i), 2. * i + 1.};
        mean << sample;
        variance << sample;
        covariance << sample;
        batch << double(i);
        autocorr << double(i);
    }
    auto means = mean.result();
    means.reduce(reducer);
    require(means.valid() == setup.have_result, "mean result survives on wrong rank");
    require(rank != 1 || (means.count() == 10 && means.mean()(0) == 5.5 && means.mean()(1) == 12.),
            "unequal runs lost counts or means");
    auto vars = variance.result();
    vars.reduce(reducer);
    require(vars.valid() == setup.have_result, "variance result survives on wrong rank");
    require(rank != 1 || (vars.count() == 10 && std::abs(vars.var()(0) - 55./6.) < 1e-12
                         && std::abs(vars.stderror()(0) - std::sqrt(11./12.)) < 1e-12),
            "reduced variance disagrees with raw samples");
    auto cov = covariance.result();
    cov.reduce(reducer);
    require(cov.valid() == setup.have_result, "covariance result survives on wrong rank");
    require(rank != 1 || (cov.count() == 10 && std::abs(cov.cov()(0,1) - 55./3.) < 1e-12),
            "reduced covariance disagrees with raw samples");
    auto batches = batch.result();
    batches.reduce(reducer);
    require(batches.valid() == setup.have_result, "batch result survives on wrong rank");
    require(rank != 1 || (batches.count() == 10 && batches.mean()(0) == 5.5 && batches.count2() == 18.),
            "reduction mixed independent partial batches");
    auto levels = autocorr.result();
    auto common = -reducer.get_max(-int64_t(levels.nlevel()));
    levels.reduce(reducer);
    require(levels.valid() == setup.have_result, "autocorrelation result survives on wrong rank");
    require(rank != 1 || (levels.count() == 10 && levels.mean()(0) == 5.5
                         && levels.nlevel() == size_t(common)),
            "unequal runs invented autocorrelation levels or lost samples");
    empty_autocorrelation(rank, reducer);
    elliptic_signed_ratio(rank, reducer);
}
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank, size, status = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    try {
        require(size == 2, "run this contract with two MPI ranks");
        contract(rank);
    } catch (std::exception const& error) {
        if (!rank) std::cerr << error.what() << '\n';
        status = 1;
    }
    MPI_Finalize();
    return status;
}
