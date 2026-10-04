// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#include <alps/alea.hpp>
#include <alps/alea/mpi.hpp>
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
