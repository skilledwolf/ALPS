// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
#include <alps/alea.hpp>
#include <alps/alea/hdf5.hpp>
#include <stdexcept>
#include <vector>

template<class T> bool close(T a, T b) {
    if (a==b) return true;
    if constexpr (std::is_floating_point_v<T>) if (std::isnan(a) && std::isnan(b)) return true;
    return std::abs(a-b) <= 1e-12 * std::max({1.,std::abs(a),std::abs(b)});
}
bool close(alps::alea::complex_op<double> a, alps::alea::complex_op<double> b) {
    return close(a.rere(),b.rere()) && close(a.reim(),b.reim())
        && close(a.imre(),b.imre()) && close(a.imim(),b.imim());
}
template<class A, class B> bool close_array(A const& a, B const& b) {
    if (a.rows()!=b.rows() || a.cols()!=b.cols()) return false;
    for (Eigen::Index i=0;i<a.rows();++i) for (Eigen::Index j=0;j<a.cols();++j)
        if (!close(a(i,j),b(i,j))) return false;
    return true;
}
template<class R> bool equivalent(R const& a, R const& b) {
    using tr=alps::alea::traits<R>;
    if (a.count()!=b.count() || !close_array(a.mean(),b.mean())) return false;
    if constexpr (tr::HAVE_VAR)
        if (a.count2()!=b.count2() || !close_array(a.var(),b.var())) return false;
    if constexpr (tr::HAVE_COV)
        if (!close_array(a.cov(),b.cov())) return false;
    if constexpr (tr::HAVE_BATCH)
        if (a.store().batch()!=b.store().batch() || a.store().count()!=b.store().count()) return false;
    if constexpr (tr::HAVE_TAU) {
        if (a.nlevel()!=b.nlevel()) return false;
        for (size_t i=0;i<a.nlevel();++i) if (!equivalent(a.level(i),b.level(i))) return false;
    }
    return true;
}

template<class Accumulator>
void check(alps::alea::hdf5_serializer& reader, char const* name, Accumulator const& expected) {
    decltype(expected.result()) actual;
    deserialize(reader, name, actual);
    // Released raw-moment arithmetic and centered accumulation round differently.
    // Counts and retained bins must remain exact; inferred moments must agree.
    if (!equivalent(actual,expected.result())) throw std::runtime_error(name);
}

int main(int argc, char** argv) {
    if (argc != 2) return 1;
    using namespace alps::alea;
    mean_acc<double> mean(3);
    var_acc<double> variance(3, 3);
    cov_acc<double> covariance(3, 3);
    autocorr_acc<double> autocorr(3, 3);
    batch_acc<double> batch(3, 16, 3);
    cov_acc<std::complex<double>, circular_var> circular(2, 3);
    cov_acc<std::complex<double>, elliptic_var> elliptic(2, 3);
    for (int i = 0; i < 67; ++i) {
        std::vector<double> values{double(i), double(i), double(i % 7)};
        auto data = make_adapter(values);
        mean << data; variance << data; covariance << data; autocorr << data; batch << data;
        std::vector<std::complex<double>> z{{double(i), double(i % 7)}, {double(i), -double(i % 7)}};
        circular << make_adapter(z); elliptic << make_adapter(z);
    }
    alps::hdf5::archive archive(argv[1], "r");
    hdf5_serializer reader(archive, "/results");
    check(reader, "mean", mean);
    check(reader, "variance", variance);
    check(reader, "covariance", covariance);
    check(reader, "autocorr", autocorr);
    check(reader, "batch", batch);
    check(reader, "circular", circular);
    check(reader, "elliptic", elliptic);
    archive.close();
}
