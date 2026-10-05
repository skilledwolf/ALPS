// SPDX-License-Identifier: MIT
#include <alps/alea.hpp>
#include <alps/alea/checkpoint.hpp>
#include <alps/alea/hdf5.hpp>
#include <alps/hdf5/vector.hpp>
#include <boost/filesystem.hpp>
#include <cmath>
#include <complex>
#include <iostream>
#include <stdexcept>
#include <type_traits>
#include <vector>
#include <array>
#include <functional>
#include <variant>

namespace {
namespace aa = alps::alea;
void require(bool condition, char const* message) {
    if (!condition) throw std::runtime_error(message);
}
template<class F> void rejects(F operation) {
    bool failed = false;
    try { operation(); } catch (std::exception const&) { failed = true; }
    require(failed, "invalid statistics operation was accepted");
}
template<class E, class F> void throws_estimate(F operation) {
    try { operation(); } catch (E const&) { return; }
    throw std::runtime_error("statistics accessor did not report its estimate limitation");
}

void mean_tests() {
    // Independent analytic oracles: scalar pooled t^2 and the inverse of a
    // two-dimensional covariance. Neither uses ALEA's diagonalization.
    aa::cov_result<double> left(aa::cov_data<double>(2)), right(aa::cov_data<double>(2));
    left.store().count()=5; left.store().count2()=5;
    left.store().data() << 2,1; left.store().data2() << 4,1,1,2;
    right.store().count()=8; right.store().count2()=8;
    right.store().data() << -1,3; right.store().data2() << 3,-.5,-.5,5;
    aa::column<double> zero(2); zero.setZero();
    auto one = aa::test_mean(left, zero);
    require(std::abs(one.score()-15./7) < 1e-12 &&
            one.dist().degrees_of_freedom1()==2 && one.dist().degrees_of_freedom2()==3,
            "one-sample covariance test disagrees with inverse-matrix oracle");
    auto two = aa::test_mean(left, right);
    double expected = 10./22 * 40./13 * 5951./1590.75;
    require(std::abs(two.score()-expected) < 1e-12 &&
            two.dist().degrees_of_freedom1()==2 && two.dist().degrees_of_freedom2()==10,
            "two-sample covariance pooling or degrees of freedom are incorrect");
    Eigen::Matrix2d rotation; rotation << .6,-.8,.8,.6;
    auto rotated = left;
    rotated.store().data() = rotation*left.mean();
    rotated.store().data2() = rotation*left.cov()*rotation.transpose();
    require(std::abs(aa::test_mean(rotated,zero).score()-one.score()) < 1e-12,
            "mean test is not invariant under orthogonal basis changes");
    auto mixed_units = left;
    Eigen::Matrix2d units = Eigen::Matrix2d::Zero(); units(0,0)=1e5; units(1,1)=1e-5;
    mixed_units.store().data() = units*left.mean();
    mixed_units.store().data2() = units*left.cov()*units;
    require(std::abs(aa::test_mean(mixed_units,zero).score()-one.score()) < 1e-12,
            "different component units changed covariance rank or significance");
    aa::column<double> diff_units(2), variance_units(2);
    diff_units << 1e5,1e-5; variance_units << 1e10,1e-10;
    require(std::abs(aa::t2_test(diff_units,variance_units,10).score()-.9) < 1e-12,
            "different component units changed a diagonal mean test");
    aa::var_result<double> a(aa::var_data<double>(1)), b(aa::var_data<double>(1));
    a.store().count()=5; a.store().count2()=5; a.store().data()(0)=2; a.store().data2()(0)=4;
    b.store().count()=8; b.store().count2()=8; b.store().data()(0)=-1; b.store().data2()(0)=3;
    auto scalar = aa::test_mean(a,b);
    require(std::abs(scalar.score()-3960./481) < 1e-12 &&
            scalar.dist().degrees_of_freedom2()==11,
            "two-sample scalar test disagrees with pooled Student t squared");
    // Fractional effective sample counts must remain doubles.
    a.store().count2()=10;
    aa::column<double> scalar_zero(1); scalar_zero.setZero();
    auto weighted = aa::test_mean(a,scalar_zero);
    require(std::abs(weighted.score()-2.5) < 1e-12 &&
            weighted.dist().degrees_of_freedom2()==1.5,
            "weighted mean test truncated the effective sample count");
    auto scaled = a;
    scaled.store().data() *= 1e-8;
    scaled.store().data2() *= 1e-16;
    require(std::abs(aa::test_mean(scaled,scalar_zero).score()-weighted.score()) < 1e-12,
            "small nonzero variance was mistaken for a deterministic result");
    aa::column<double> huge(1), huge_variance(1);
    huge(0)=1e155; huge_variance(0)=1e308;
    require(std::abs(aa::t2_test(huge,huge_variance,10).score()-100) < 1e-12,
            "finite standardized mean test overflowed while squaring its difference");
    // Noncircular complex data must match a joint real/imaginary test.
    aa::batch_acc<std::complex<double>> complex(1,8,1);
    aa::batch_acc<double> real(2,8,1);
    for (size_t i=0; i<7; ++i) {
        double x=double(i)-2, y=double((i*i)%5)+.2*x;
        complex << std::vector<std::complex<double>>{{x,y}};
        real << std::vector<double>{x,y};
    }
    aa::column<std::complex<double>> complex_zero(1); complex_zero.setZero();
    auto complex_test = aa::test_mean(complex.result(),complex_zero);
    auto real_test = aa::test_mean(real.result(),zero);
    require(std::abs(complex_test.score()-real_test.score()) < 1e-12 &&
            complex_test.dist().degrees_of_freedom1()==2,
            "complex mean test lost real/imaginary covariance or dimension");
    aa::cov_acc<std::complex<double>,aa::elliptic_var> elliptic(1);
    aa::var_acc<std::complex<double>,aa::elliptic_var> elliptic_var(1);
    aa::cov_acc<std::complex<double>> circular(1);
    for (size_t i=0; i<7; ++i) {
        auto value=std::vector<std::complex<double>>{{double(i)-2,double((i*i)%5)+.2*(double(i)-2)}};
        elliptic << value; elliptic_var << value; circular << value;
    }
    require(std::abs(aa::test_mean(elliptic.result(),complex_zero).score()-real_test.score()) < 1e-12 &&
            std::abs(aa::test_mean(elliptic_var.result(),complex_zero).score()-real_test.score()) < 1e-12,
            "elliptic result mean test differs from joint real/imaginary oracle");
    rejects([&] { aa::test_mean(circular.result(),complex_zero); });
    aa::var_acc<double> missing_covariance(2);
    for (size_t i=0; i<7; ++i) missing_covariance << std::vector<double>{double(i),double(i*i)};
    rejects([&] { aa::test_mean(missing_covariance.result(),zero); });
    a.store().data2()(0)=0;
    require(aa::test_mean(a,a.mean()).pvalue()==1 && aa::test_mean(a,scalar_zero).pvalue()==0,
            "deterministic equality or mismatch produced an undefined probability");
    a.store().data2()(0)=-1;
    rejects([&] { aa::test_mean(a,scalar_zero); });
    a.store().data2()(0)=NAN;
    rejects([&] { aa::test_mean(a,scalar_zero); });
    rejects([&] { aa::test_mean(aa::var_result<double>(),scalar_zero); });
    a.store().count2()=.25;
    rejects([&] { aa::test_mean(a,scalar_zero); });
    rejects([&] { aa::t2_test(aa::column<double>{},aa::column<double>{},10); });
    auto asymmetric = left; asymmetric.store().data2()(0,1)=.5;
    rejects([&] { aa::test_mean(asymmetric,zero); });
    asymmetric = mixed_units; asymmetric.store().data2()(0,1)*=2;
    rejects([&] { aa::test_mean(asymmetric,zero); });
    rejects([&] { aa::test_mean(left,scalar_zero); });
    rejects([&] { aa::t2_test(zero,aa::column<double>::Ones(2),.5); });
    aa::cov_acc<double> undersampled(2);
    undersampled << std::vector<double>{1,0} << std::vector<double>{0,1};
    rejects([&] { aa::test_mean(undersampled.result(),zero); });
    aa::fisher_f_distribution f(1,1);
    double lower = 2/std::acos(-1.)*std::atan(std::sqrt(3.));
    require(std::abs(f.cdf(3)-lower) < 1e-14 && std::abs(f.ccdf(3)-(1-lower)) < 1e-14 &&
            f.ccdf(INFINITY)==0 && f.cdf(INFINITY)==1 &&
            std::isnan(aa::fisher_f_distribution(0,1).cdf(1)),
            "F distribution tails disagree with the analytic F(1,1) law");
    require(aa::t2_result(0,4,10).pvalue()==1,
            "standard mean-test p-value incorrectly counts lower-tail variance anomalies");
    auto tiny_tail = aa::fisher_f_distribution(2,2).ccdf(1e308);
    require(tiny_tail > 0 && std::abs(tiny_tail/1e-308-1) < 1e-13,
            "F distribution overflow discarded a representable extreme upper tail");
}

// A two-participant sum reducer that deliberately defers writes until commit.
// The peer's primitive contributions are recorded separately from the target.
struct test_reducer : aa::reducer {
    using buffer = std::variant<std::vector<double>, std::vector<int32_t>,
                                std::vector<int64_t>, std::vector<uint64_t>>;
    aa::reducer_setup setup;
    std::vector<int64_t> maxima;
    std::vector<buffer> const* peer;
    bool fail_commit;
    mutable size_t query = 0, field = 0;
    mutable std::vector<buffer> recorded;
    mutable std::vector<std::function<void()>> pending;
    test_reducer(aa::reducer_setup setup, std::vector<int64_t> maxima,
                 std::vector<buffer> const* peer = nullptr, bool fail = false)
        : setup(setup), maxima(std::move(maxima)), peer(peer), fail_commit(fail) {}
    aa::reducer_setup get_setup() const override { return setup; }
    int64_t get_max(int64_t value) const override {
        require(query < maxima.size() && value <= maxima[query], "invalid reduction consensus query");
        return maxima[query++];
    }
    template<class T> void sum(aa::view<T> data) const {
        if (!peer) {
            std::vector<T> copy(data.size());
            for (size_t i=0; i<data.size(); ++i) copy[i] = data.data()[i];
            recorded.emplace_back(std::move(copy));
        } else {
            require(field < peer->size(), "extra reduction payload");
            auto const& contribution = std::get<std::vector<T>>((*peer)[field]);
            require(contribution.size() == data.size(), "inconsistent reduction payload shapes");
            pending.emplace_back([data, contribution]() mutable {
                for (size_t i=0; i<data.size(); ++i) data.data()[i] += contribution[i];
            });
        }
        ++field;
    }
    void reduce(aa::view<double> data) const override { sum(data); }
    void reduce(aa::view<int32_t> data) const override { sum(data); }
    void reduce(aa::view<int64_t> data) const override { sum(data); }
    void reduce(aa::view<uint64_t> data) const override { sum(data); }
    void commit() const override {
        for (auto const& operation : pending) operation();
        if (fail_commit) throw std::runtime_error("deferred reduction failed");
    }
};
std::vector<int64_t> dimensions(size_t size) { return {0,0,int64_t(size),-int64_t(size)}; }
template<class R> R reduced(R const& left, R const& right, std::vector<int64_t> maxima) {
    test_reducer capture({1,2,false}, maxima);
    auto peer = right;
    peer.reduce(capture);
    require(!peer.valid(), "nonrecipient retained a reduction result");
    test_reducer sum({0,2,true}, maxima, &capture.recorded);
    auto result = left;
    result.reduce(sum);
    require(capture.query == maxima.size() && sum.query == maxima.size(),
            "missing reduction consensus query");
    require(sum.field == capture.recorded.size(), "missing reduction payload");
    return result;
}
template<class R> void reduction_failure(R const& left, R const& right,
                                         std::vector<int64_t> maxima) {
    test_reducer capture({1,2,false}, maxima);
    auto peer = right;
    peer.reduce(capture);
    test_reducer sum({0,2,true}, maxima, &capture.recorded, true);
    auto result = left;
    rejects([&] { result.reduce(sum); });
    require(result == left, "failed deferred reduction mutated the original result");
}
template<class T> std::vector<T> reduction_sample(size_t i, int run) {
    double a = double(int(i%11)-5) + .25*double(i%3) + 3*run;
    double b = double(int((7*i)%13)-6) - .5*double(i%4) - 2*run;
    if constexpr (std::is_same_v<T, std::complex<double>>)
        return {T(a,.5*b), T(b,-.25*a+run)};
    else return {a,b};
}
template<class T> struct bin { uint64_t count=0; std::array<T,2> sum{}; };
template<class T> std::vector<bin<T>> raw_bins(int run, size_t count, size_t width) {
    std::vector<bin<T>> bins;
    for (size_t i=0; i<count; ++i) {
        if (bins.empty() || bins.back().count == width) bins.emplace_back();
        ++bins.back().count;
        auto value = reduction_sample<T>(i,run);
        for (size_t j=0; j<2; ++j) bins.back().sum[j] += value[j];
    }
    return bins;
}
template<class T> T conjugate(T value) {
    if constexpr (std::is_same_v<T, std::complex<double>>) return std::conj(value);
    else return value;
}
template<class T> struct bin_oracle {
    uint64_t count=0;
    double count2=0;
    std::array<T,2> mean{};
    std::array<std::array<T,2>,2> cov{};
    explicit bin_oracle(std::vector<bin<T>> const& bins) {
        for (auto const& b : bins) {
            count += b.count;
            count2 += double(b.count)*b.count;
            for (size_t j=0; j<2; ++j) mean[j] += b.sum[j];
        }
        for (auto& m : mean) m /= double(count);
        // Independent two-pass weighted covariance of the original bin means.
        for (auto const& b : bins) if (b.count) {
            for (size_t j=0; j<2; ++j) for (size_t k=0; k<2; ++k)
                cov[j][k] += double(b.count) * (b.sum[j]/double(b.count)-mean[j])
                            * conjugate(b.sum[k]/double(b.count)-mean[k]);
        }
        for (auto& row : cov) for (auto& c : row) c /= double(count)-count2/count;
    }
    template<class R> void check_variance(R const& result) const {
        require(result.count() == count && result.count2() == count2,
                "independent reduction lost samples or original bin weights");
        for (size_t j=0; j<2; ++j) {
            require(std::abs(result.mean()(j)-mean[j]) < 1e-11, "independent reduction changed mean");
            require(std::abs(result.var()(j)-std::real(cov[j][j])) < 1e-10,
                    "independent bin variance disagrees with two-pass oracle");
            auto error = std::sqrt(std::real(cov[j][j])*count2/(double(count)*count));
            require(std::abs(result.stderror()(j)-error) < 1e-11,
                    "independent bin error disagrees with two-pass oracle");
        }
    }
};
template<class T> void independent_reductions() {
    aa::mean_acc<T> mean_left(2), mean_right(2);
    aa::var_acc<T> var_left(2,3), var_right(2,3);
    aa::cov_acc<T> cov_left(2,3), cov_right(2,3);
    aa::autocorr_acc<T> auto_left(2,3,2), auto_right(2,3,2);
    aa::batch_acc<T> batch_left(2,4,3), batch_right(2,8,5);
    for (size_t i=0; i<53; ++i) {
        auto right = reduction_sample<T>(i,1);
        mean_right << right; var_right << right; cov_right << right;
        auto_right << right; batch_right << right;
        if (i<19) {
            auto left = reduction_sample<T>(i,0);
            mean_left << left; var_left << left; cov_left << left;
            auto_left << left; batch_left << left;
        }
    }
    auto bins = raw_bins<T>(0,19,3), peer_bins = raw_bins<T>(1,53,3);
    bins.insert(bins.end(), peer_bins.begin(), peer_bins.end());
    bin_oracle<T> oracle(bins);
    auto mean = reduced(mean_left.result(), mean_right.result(), dimensions(2));
    reduction_failure(mean_left.result(), mean_right.result(), dimensions(2));
    require(mean.count() == 72, "unequal-run mean lost samples");
    for (size_t j=0; j<2; ++j) require(std::abs(mean.mean()(j)-oracle.mean[j]) < 1e-11,
                                    "unequal-run mean ignored run lengths");
    auto maxima = dimensions(2); maxima.insert(maxima.end(), {0,0});
    oracle.check_variance(reduced(var_left.result(), var_right.result(), maxima));
    reduction_failure(var_left.result(), var_right.result(), maxima);
    auto covariance = reduced(cov_left.result(), cov_right.result(), maxima);
    reduction_failure(cov_left.result(), cov_right.result(), maxima);
    oracle.check_variance(covariance);
    for (size_t j=0; j<2; ++j) for (size_t k=0; k<2; ++k)
        require(std::abs(covariance.cov()(j,k)-oracle.cov[j][k]) < 1e-10,
                "independent reduction changed cross-covariance");
    auto left = batch_left.result(), right = batch_right.result();
    std::vector<bin<T>> original_bins;
    for (auto const* result : {&left,&right}) for (size_t i=0; i<result->num_batches(); ++i) {
        bin<T> b; b.count = result->store().count()(i);
        for (size_t j=0; j<2; ++j) b.sum[j] = result->store().batch()(j,i);
        original_bins.push_back(b);
    }
    maxima = dimensions(2); maxima.insert(maxima.end(), {0,0,8});
    auto batch = reduced(left, right, maxima);
    reduction_failure(left, right, maxima);
    require(batch.count() == 72, "partial independent batches lost samples");
    for (size_t j=0; j<2; ++j) require(std::abs(batch.mean()(j)-oracle.mean[j]) < 1e-11,
                                     "partial independent batches changed raw mean");
    bin_oracle<T> batch_oracle(original_bins);
    batch_oracle.check_variance(batch);
    for (size_t j=0; j<2; ++j) for (size_t k=0; k<2; ++k)
        require(std::abs(batch.cov()(j,k)-batch_oracle.cov[j][k]) < 1e-10,
                "reduction mixed independent run bins");
    auto al = auto_left.result(), ar = auto_right.result();
    auto common = std::min(al.nlevel(),ar.nlevel());
    require(al.nlevel() != ar.nlevel(), "autocorrelation fixture needs unequal hierarchy depths");
    maxima = dimensions(2); maxima.insert(maxima.end(), {0,0,-int64_t(common)});
    maxima.insert(maxima.end(), 2*common, 0);
    auto autocorr = reduced(al, ar, maxima);
    reduction_failure(al, ar, maxima);
    require(autocorr.nlevel() == common, "reduction invented unavailable coarse levels");
    for (size_t level=0, width=3; level<common; ++level, width*=2) {
        bins = raw_bins<T>(0,19,width); peer_bins = raw_bins<T>(1,53,width);
        bins.insert(bins.end(), peer_bins.begin(), peer_bins.end());
        bin_oracle<T>(bins).check_variance(autocorr.level(level));
    }
}
template<class T> void empty_autocorrelation_reduction() {
    aa::autocorr_acc<T> populated(2), empty(2,9,3);
    constexpr size_t samples = 8195;
    for (size_t i=0; i<samples; ++i) populated << reduction_sample<T>(i/16,0);
    auto original = populated.result();
    auto selected = original.find_level(1024);
    require(selected > 0 && original.tau_available(), "empty-run fixture needs coarse error estimates");
    auto raw = [](size_t width) {
        std::vector<bin<T>> bins;
        for (size_t i=0; i<samples; ++i) {
            if (bins.empty() || bins.back().count == width) bins.emplace_back();
            ++bins.back().count;
            auto value = reduction_sample<T>(i/16,0);
            for (size_t j=0; j<2; ++j) bins.back().sum[j] += value[j];
        }
        return bin_oracle<T>(bins);
    };
    auto base = raw(1), coarse = raw(size_t(1)<<selected);
    auto check = [&](aa::autocorr_result<T> const& result) {
        require(result.count() == samples && result.nlevel() == original.nlevel()
             && result.count2() == original.count2()
             && result.observations() == original.observations(),
                "empty run changed autocorrelation depth or effective sample count");
        base.check_variance(result.level(0));
        coarse.check_variance(result.level(selected));
        for (size_t level=0; level<original.nlevel(); ++level)
            require(result.level(level).count() == original.level(level).count()
                 && result.level(level).count2() == original.level(level).count2(),
                    "empty run changed autocorrelation level weights");
        for (size_t j=0; j<2; ++j) {
            auto variance = std::real(coarse.cov[j][j]);
            auto error = std::sqrt(variance*coarse.count2/(double(samples)*samples));
            auto tau = .5*(coarse.count2/samples)*variance/std::real(base.cov[j][j])-.5;
            require(std::abs(result.mean()(j)-coarse.mean[j]) < 1e-11
                 && std::abs(result.var()(j)-variance) < 1e-10
                 && std::abs(result.stderror()(j)-error) < 1e-11
                 && std::abs(result.tau()(j)-tau) < 1e-10,
                    "empty-run autocorrelation statistics disagree with raw-bin oracle");
        }
    };
    check(original);
    auto maxima = dimensions(2);
    maxima.insert(maxima.end(), {0,0,-int64_t(original.nlevel())});
    maxima.insert(maxima.end(), 2*original.nlevel(), 0);
    check(reduced(original, empty.result(), maxima));
    check(reduced(empty.result(), original, maxima));
    reduction_failure(original, empty.result(), maxima);

    aa::autocorr_result<T> multiple_empty(3);
    for (size_t i=0; i<multiple_empty.nlevel(); ++i)
        multiple_empty.level(i) = aa::var_result<T>(aa::var_data<T>(2));
    maxima = dimensions(2);
    maxima.insert(maxima.end(), {0,0,std::numeric_limits<int64_t>::min(),0,0});
    auto all_empty = reduced(empty.result(), multiple_empty, maxima);
    require(all_empty.valid() && all_empty.size() == 2 && all_empty.count() == 0
         && all_empty.nlevel() == 1 && all_empty.level(0).count2() == 0.,
            "all-empty autocorrelation reduction invented samples or retained redundant levels");
    multiple_empty.level(2).store().count2() = 1.;
    maxima = dimensions(2); maxima.push_back(1);
    test_reducer malformed({0,1,true},maxima);
    rejects([&] { multiple_empty.reduce(malformed); });
    require(malformed.field == 0 && multiple_empty.nlevel() == 3
         && multiple_empty.level(2).count2() == 1.,
            "discarded malformed empty level was hidden or changed");
}
void reduction_limits() {
    uint64_t a = uint64_t(1)<<62, b = a+1, total = a+b;
    aa::mean_data<double> ml(1), mr(1);
    ml.count() = a; mr.count() = b; ml.data()(0) = mr.data()(0) = 1.;
    auto mean = reduced(aa::mean_result<double>(ml), aa::mean_result<double>(mr), dimensions(1));
    require(mean.count() == total && mean.mean()(0) == 1., "reduction truncated unsigned sample count");
    aa::var_data<double> vl(1), vr(1);
    vl.count() = a; vr.count() = b; vl.count2() = double(a); vr.count2() = double(b);
    vl.data()(0) = vr.data()(0) = 1.;
    auto maxima = dimensions(1); maxima.insert(maxima.end(), {0,0});
    auto variance = reduced(aa::var_result<double>(vl), aa::var_result<double>(vr), maxima);
    require(variance.count() == total && variance.mean()(0) == 1. && variance.var()(0) == 0.,
            "variance reduction truncated unsigned sample count");
    aa::cov_data<double> cl(1), cr(1);
    cl.count() = a; cr.count() = b; cl.count2() = double(a); cr.count2() = double(b);
    cl.data()(0) = cr.data()(0) = 1.;
    auto covariance = reduced(aa::cov_result<double>(cl), aa::cov_result<double>(cr), maxima);
    require(covariance.count() == total && covariance.mean()(0) == 1. && covariance.cov()(0,0) == 0.,
            "covariance reduction truncated unsigned sample count");

    aa::batch_acc<double> empty(2,4,3), populated(2,8,3);
    for (size_t i=0; i<19; ++i) populated << reduction_sample<double>(i,0);
    auto original = populated.result();
    maxima = dimensions(2); maxima.insert(maxima.end(), {0,0,8});
    auto combined = reduced(empty.result(), original, maxima);
    require(combined.count() == original.count() && combined.count2() == original.count2()
         && combined.mean() == original.mean() && combined.cov().isApprox(original.cov(),1e-12),
            "empty independent run changed retained batch statistics");

    aa::batch_data<double> bad_bins(2,2);
    bad_bins.batch()(0,0) = 1.;
    aa::batch_result<double> malformed(bad_bins), previous = malformed;
    maxima = dimensions(2); maxima.push_back(1);
    test_reducer invalid_bins({0,1,true},maxima);
    rejects([&] { malformed.reduce(invalid_bins); });
    require(malformed == previous && invalid_bins.field == 0,
            "zero-weight nonzero bin was reduced or mutated");

    aa::cov_acc<double> acc(2,3);
    for (size_t i=0; i<19; ++i) acc << reduction_sample<double>(i,0);
    auto target = acc.result(), before = target;
    test_reducer wrong_size({0,2,true},{0,0,3,-2});
    rejects([&] { target.reduce(wrong_size); });
    require(target == before && wrong_size.field == 0,
            "shape disagreement exposed data or changed result");
    target.store().count2() = 0.;
    before = target;
    maxima = dimensions(2); maxima.push_back(1);
    test_reducer invalid_weights({0,1,true},maxima);
    rejects([&] { target.reduce(invalid_weights); });
    require(target == before && invalid_weights.field == 0,
            "invalid squared-weight count reached reduction or mutated result");
    test_reducer invalid_setup({1,1,true},{});
    rejects([&] { target.reduce(invalid_setup); });
    require(target == before, "invalid reducer setup mutated result");
}
void elliptic_reduction() {
    using T = std::complex<double>;
    aa::var_acc<T,aa::elliptic_var> vl(2,3), vr(2,3);
    aa::cov_acc<T,aa::elliptic_var> cl(2,3), cr(2,3);
    for (size_t i=0; i<53; ++i) {
        auto right = reduction_sample<T>(i,1); vr << right; cr << right;
        if (i<19) { auto left = reduction_sample<T>(i,0); vl << left; cl << left; }
    }
    auto bins = raw_bins<T>(0,19,3), peer = raw_bins<T>(1,53,3);
    bins.insert(bins.end(),peer.begin(),peer.end());
    bin_oracle<T> oracle(bins);
    auto maxima = dimensions(2); maxima.insert(maxima.end(), {0,0});
    auto variance = reduced(vl.result(),vr.result(),maxima);
    auto covariance = reduced(cl.result(),cr.result(),maxima);
    require(variance.count() == oracle.count && covariance.count() == oracle.count
         && variance.count2() == oracle.count2 && covariance.count2() == oracle.count2,
            "elliptic reduction lost weights");
    for (size_t j=0; j<2; ++j) for (size_t k=0; k<2; ++k) {
        double expected[4]{};
        for (auto const& b : bins) {
            auto x = b.sum[j]/double(b.count)-oracle.mean[j];
            auto y = b.sum[k]/double(b.count)-oracle.mean[k];
            expected[0] += b.count*x.real()*y.real();
            expected[1] += b.count*x.real()*y.imag();
            expected[2] += b.count*x.imag()*y.real();
            expected[3] += b.count*x.imag()*y.imag();
        }
        auto check = [&](aa::complex_op<double> value, double scale) {
            double actual[]{value.rere(),value.reim(),value.imre(),value.imim()};
            for (size_t axis=0; axis<4; ++axis)
                require(std::abs(actual[axis]-scale*expected[axis]
                       /(oracle.count-oracle.count2/oracle.count)) < 1e-10,
                        "elliptic reduction lost real/imaginary covariance");
        };
        check(covariance.cov()(j,k),1.);
        if (j==k) {
            check(variance.var()(j),1.);
            auto error = variance.stderror()(j);
            check(dot(error,error),oracle.count2/(double(oracle.count)*oracle.count));
        }
    }
}
template<class R> R roundtrip(std::string const& filename, std::string const& name, R const& result) {
    {
        alps::hdf5::archive ar(filename, "a");
        aa::hdf5_serializer codec(ar, "/results");
        serialize(codec, name, result);
        ar.close();
    }
    R restored;
    {
        alps::hdf5::archive ar(filename, "r");
        aa::hdf5_serializer codec(ar, "/results");
        deserialize(codec, name, restored);
    }
    require(result == restored, "scientific result failed exact round trip");
    return restored;
}
template<class T> std::vector<T> sample(size_t i, size_t size) {
    std::vector<T> value(size);
    for (size_t j=0; j<size; ++j) {
        auto real = double(int(i%17)-8) + .125*double(i%3) + double(j);
        if constexpr (std::is_same_v<T, std::complex<double>>)
            value[j] = T(real, .5*real - double(j));
        else value[j] = real;
    }
    return value;
}
template<class T> void resume(std::string const& filename, size_t components, size_t boundary, uint64_t base) {
    aa::batch_acc<T> uninterrupted(components, 8, base);
    for (size_t i=0; i<boundary; ++i) uninterrupted << sample<T>(i, components);
    {
        alps::hdf5::archive ar(filename, "w");
        aa::hdf5_serializer codec(ar, "/");
        serialize(codec, "state", uninterrupted);
        require(ar.extent("/state/batch/sum") == std::vector<size_t>({8, components}), "batch axes changed");
        require(ar.list_children("/state/cursor").size() == 2, "checkpoint stored redundant cursor fields");
    }
    aa::batch_acc<T> restored(1, 2, 1);
    {
        alps::hdf5::archive ar(filename, "r");
        aa::hdf5_serializer codec(ar, "/");
        deserialize(codec, "state", restored);
    }
    for (size_t i=boundary; i<281; ++i) {
        auto value = sample<T>(i, components);
        uninterrupted << value;
        restored << value;
    }
    require(uninterrupted.result() == restored.result(), "partial checkpoint changed subsequent batches");
    require(uninterrupted.result().store().count() == restored.result().store().count(),
            "partial checkpoint changed individual batch weights");
    require(uninterrupted.offset() == restored.offset(), "partial checkpoint changed time ordering");
    require(uninterrupted.cursor().current() == restored.cursor().current()
         && uninterrupted.cursor().cycle() == restored.cursor().cycle()
         && uninterrupted.cursor().level() == restored.cursor().level(), "merge cursor changed after resume");
    roundtrip(filename, "batch", restored.result());
}
void results(std::string const& filename) {
    aa::mean_acc<double> mean(1);
    aa::var_acc<double> variance(3);
    aa::cov_acc<double> covariance(3);
    aa::autocorr_acc<double> autocorr(3);
    aa::cov_acc<std::complex<double>> complex_cov(2);
    aa::cov_acc<std::complex<double>, aa::elliptic_var> elliptic_cov(2);
    for (size_t i=0; i<97; ++i) {
        mean << sample<double>(i, 1);
        variance << sample<double>(i, 3);
        covariance << sample<double>(i, 3);
        autocorr << sample<double>(i, 3);
        complex_cov << sample<std::complex<double>>(i, 2);
        elliptic_cov << sample<std::complex<double>>(i, 2);
    }
    roundtrip(filename, "mean", mean.result());
    roundtrip(filename, "variance", variance.result());
    roundtrip(filename, "covariance", covariance.result());
    roundtrip(filename, "autocorr", autocorr.result());
    roundtrip(filename, "complex-covariance", complex_cov.result());
    roundtrip(filename, "elliptic-covariance", elliptic_cov.result());
    {
        alps::hdf5::archive ar(filename, "a");
        aa::hdf5_serializer codec(ar, "/results");
        aa::batch_result<double> empty(aa::batch_data<double>(1,0));
        serialize(codec, "empty-batch", empty);
        aa::batch_result<double> restored;
        deserialize(codec, "empty-batch", restored);
        require(restored.count() == 0 && restored.size() == 1 && restored.store().num_batches() == 0,
                "empty analysis result acquired invented samples");
    }
    alps::hdf5::archive ar(filename, "r");
    require(ar.extent("/results/covariance/cov") == std::vector<size_t>({3,3}), "covariance axes changed");
    require(ar.extent("/results/elliptic-covariance/cov") == std::vector<size_t>({2,2,2,2}), "elliptic covariance lost real operator axes");
    require(ar.is_datatype<std::complex<double>>("/results/complex-covariance/cov"), "complex covariance is not canonical");
}
void unequal_merge(std::string const& filename) {
    aa::var_acc<double> left, right;
    for (auto v : {1., 2., 3.}) left << v;
    for (auto v : {4., 5., 6., 7., 8., 9., 10.}) right << v;
    auto right_result = right.result();
    auto restored = roundtrip(filename, "unequal-right", right_result);
    left << restored;
    auto merged = left.result();
    require(merged.count() == 10 && merged.mean()(0) == 5.5, "merge lost unequal sample weights");
    require(std::abs(merged.var()(0) - 55./6.) < 1e-12, "merged variance disagrees with raw samples");
    require(std::abs(merged.stderror()(0) - std::sqrt(11./12.)) < 1e-12, "merged error disagrees with raw samples");
}
void batch_reset_and_equality() {
    aa::batch_acc<double> large_base(1, 2, std::numeric_limits<uint64_t>::max());
    require(large_base.current_batch_size() == std::numeric_limits<uint64_t>::max(),
            "batch size rounded outside integer range");
    large_base << 1.;
    require(large_base.count() == 1, "large batch size did not preserve a sample");
    aa::batch_acc<double> changed(2,8,3);
    for (size_t slots : {size_t(2), size_t(16)}) {
        changed.set_num_batches(slots);
        aa::batch_acc<double> fresh(2,slots,3);
        for (size_t i=0; i<281; ++i) {
            auto value = sample<double>(i,2);
            changed << value;
            fresh << value;
        }
        require(changed.result() == fresh.result() && changed.offset() == fresh.offset(),
                "resizing batches retained the old merge cursor or offsets");
    }
    auto previous = changed.result();
    rejects([&] { changed.set_num_batches(3); });
    rejects([&] { changed.set_batch_size(0); });
    require(changed.result() == previous, "invalid batch setter destroyed accumulated samples");
    auto altered = previous.store();
    size_t nonempty = 0;
    while (nonempty < altered.num_batches() && !altered.count()(nonempty)) ++nonempty;
    require(nonempty < altered.num_batches(), "batch equality fixture requires a nonempty slot");
    --altered.count()(nonempty);
    ++altered.count()((nonempty+1)%altered.num_batches());
    aa::batch_result<double> redistributed(altered);
    require(redistributed.count() == previous.count() && redistributed != previous,
            "batch equality ignored individual sample weights");
    changed.finalize();
    changed.set_num_batches(4);
    require(!changed.valid(), "reset setter revived a finalized accumulator");
    changed.reset();
    aa::batch_acc<double> fresh(2,4,3);
    for (size_t i=0; i<97; ++i) {
        auto value = sample<double>(i,2);
        changed << value;
        fresh << value;
    }
    require(changed.result() == fresh.result(), "finalized resize reset retained old state");
}
void large_batch_weights() {
    aa::batch_data<double> data(1,2);
    uint64_t const weight = uint64_t(1)<<32;
    data.count() << weight, weight;
    data.batch() << double(weight), 3.*double(weight);
    aa::batch_result<double> result(data);
    require(result.count2() == 2.*double(weight)*double(weight)
         && result.batch_size() == double(weight) && result.observations() == 2.,
            "large batch weights corrupted effective sample counts");
    require(result.mean()(0) == 2. && result.var()(0) == 2.
         && result.cov()(0,0) == 2. && result.stderror()(0) == 1.,
            "large integer batch weights overflowed squared-weight counts");
}
void wrong_sized_append() {
    aa::batch_acc<double> batch(1,8,1);
    for (size_t i=0; i<8; ++i) batch << double(i);
    auto control = batch;
    rejects([&] { batch << std::vector<double>{1.,2.}; });
    require(batch.result() == control.result() && batch.offset() == control.offset()
         && batch.cursor().current() == control.cursor().current()
         && batch.cursor().level() == control.cursor().level()
         && batch.cursor().cycle() == control.cursor().cycle(),
            "wrong-sized sample changed a full batch before rejection");
    for (size_t i=8; i<37; ++i) { batch << double(i); control << double(i); }
    require(batch.result() == control.result(), "rejected batch sample changed later continuation");
    for (size_t boundary : {size_t(0), size_t(2)}) {
        aa::autocorr_acc<double> acc(1,3);
        for (size_t i=0; i<boundary; ++i) acc << double(i);
        auto unchanged = acc;
        rejects([&] { acc << std::vector<double>{1.,2.}; });
        require(acc.count() == unchanged.count() && acc.nlevel() == unchanged.nlevel(),
                "wrong-sized sample changed autocorrelation hierarchy before rejection");
        for (size_t i=boundary; i<97; ++i) { acc << double(i); unchanged << double(i); }
        require(acc.result() == unchanged.result(), "rejected autocorrelation sample changed later estimates");
    }
}
void covariance_transform(std::string const& filename) {
    aa::cov_acc<double> acc(2);
    for (int i=0; i<64; ++i) acc << std::vector<double>{double(i), double(i)};
    auto restored = roundtrip(filename, "correlated", acc.result());
    Eigen::Matrix<double,1,2> difference;
    difference << 1., -1.;
    auto transformed = aa::transform(aa::linear_prop(), aa::linear_transformer<double>(difference), restored);
    require(transformed.size() == 1 && transformed.mean()(0) == 0., "non-square transform has incorrect orientation");
    require(std::abs(transformed.var()(0)) < 1e-20, "serialization discarded x-y covariance");
}
struct signed_batch_oracle {
    aa::column<double> mean;
    Eigen::MatrixXd covariance;
    double count, count2;
};
signed_batch_oracle signed_oracle(aa::batch_data<double> const& data) {
    signed_batch_oracle expected{aa::column<double>::Zero(data.size()),
        Eigen::MatrixXd::Zero(data.size(), data.size()), 0., 0.};
    for (size_t i=0; i<data.num_batches(); ++i) {
        auto weight = double(data.count()(i));
        expected.count += weight;
        expected.count2 += weight*weight;
        expected.mean += data.batch().col(i);
    }
    expected.mean /= expected.count;
    for (size_t i=0; i<data.num_batches(); ++i) {
        auto weight = double(data.count()(i));
        if (!weight) continue;
        aa::column<double> deviation = data.batch().col(i)/weight - expected.mean;
        expected.covariance += weight * deviation * deviation.transpose();
    }
    expected.covariance /= expected.count - expected.count2/expected.count;
    return expected;
}
void elliptic_signed_ratios(std::string const& filename) {
    using complex = std::complex<double>;
    using accumulator = aa::var_acc<complex,aa::elliptic_var>;
    accumulator left(2,3), right(2,3);
    std::vector<bin<complex>> bins;
    for (int run=0; run<2; ++run) {
        for (size_t i=0; i<(run ? 53 : 19); ++i) {
            if (i%3 == 0) bins.emplace_back();
            auto sign = i%5 == 0 ? -1. : 1.;
            std::vector<complex> value{{sign*(1. + .125*(i%7) + .5*run),sign},
                                      {sign*(-2. + .0625*(i%9) - .25*run),sign}};
            (run ? right : left) << value;
            ++bins.back().count;
            for (size_t j=0; j<2; ++j) bins.back().sum[j] += value[j];
        }
    }
    auto maxima = dimensions(2); maxima.insert(maxima.end(), {0,0});
    auto pooled = reduced(left.result(),right.result(),maxima);
    pooled = roundtrip(filename,"elliptic-signed-joint",pooled);
    auto before = pooled;
    auto ratio = aa::ratio_real_imag(pooled);
    require(pooled == before && ratio.count() == pooled.count() && ratio.count2() == pooled.count2(),
            "componentwise ratio changed its joint state or bin weights");
    for (size_t j=0; j<2; ++j) {
        aa::batch_data<double> pairs(2,bins.size());
        for (size_t i=0; i<bins.size(); ++i) {
            pairs.count()(i) = bins[i].count;
            pairs.batch()(0,i) = bins[i].sum[j].real();
            pairs.batch()(1,i) = bins[i].sum[j].imag();
        }
        auto oracle = signed_oracle(pairs);
        Eigen::Vector2d gradient(1./oracle.mean(1), -oracle.mean(0)/std::pow(oracle.mean(1),2));
        auto variance = double((gradient.transpose()*oracle.covariance*gradient)(0,0));
        auto error = std::sqrt(variance*oracle.count2/(oracle.count*oracle.count));
        require(std::abs(ratio.mean()(j)-oracle.mean(0)/oracle.mean(1)) < 1e-13
             && std::abs(ratio.var()(j)-variance) < 1e-11
             && std::abs(ratio.stderror()(j)-error) < 1e-12,
                "componentwise signed ratio disagrees with independent raw covariance");
        auto diagonal_only = gradient(0)*gradient(0)*oracle.covariance(0,0)
                           + gradient(1)*gradient(1)*oracle.covariance(1,1);
        require(diagonal_only > 2*variance, "elliptic fixture does not distinguish sign covariance");
    }
    roundtrip(filename,"elliptic-signed-ratio",ratio);

    accumulator empty(2), single(1), zero_sign(1), constant(2);
    auto unavailable = aa::ratio_real_imag(empty.result());
    require(unavailable.size() == 2 && unavailable.count() == 0 && unavailable.count2() == 0.
         && unavailable.mean().array().isNaN().all() && unavailable.stderror().array().isNaN().all(),
            "empty ratio lost its shape or invented an estimate");
    single << complex(4.,-1.);
    unavailable = aa::ratio_real_imag(single.result());
    require(unavailable.count() == 1 && unavailable.mean()(0) == -4.
         && std::isnan(unavailable.stderror()(0)), "single observation acquired an independent error");
    zero_sign << complex(2.,1.) << complex(-3.,-1.);
    auto zero = zero_sign.result();
    throws_estimate<std::domain_error>([&] { aa::ratio_real_imag(zero); });
    require(zero == zero_sign.result(), "undefined ratio changed the original result");
    for (size_t i=0; i<193; ++i) {
        auto sign = i%4 == 0 ? -1. : 1.;
        constant << std::vector<complex>{{2*sign,sign},{3.,1.}};
    }
    auto exact = aa::ratio_real_imag(constant.result());
    require(exact.mean()(0) == 2. && exact.mean()(1) == 3.
         && exact.stderror().isZero(1e-7), "constant correlated ratios acquired spurious errors");

    aa::var_data<complex,aa::elliptic_var> data(1);
    data.count() = 100; data.count2() = 100.; data.data()(0) = complex(2.,1.);
    auto epsilon = std::numeric_limits<double>::epsilon();
    data.data2()(0) = aa::complex_op<double>(4.-4*epsilon,2.,2.,1.);
    auto clamped = aa::ratio_real_imag(aa::var_result<complex,aa::elliptic_var>(data));
    require(clamped.var()(0) == 0., "cancelled signed variance was not bounded by roundoff");
    data.data2()(0).rere() = 3.;
    throws_estimate<std::domain_error>([&] {
        aa::ratio_real_imag(aa::var_result<complex,aa::elliptic_var>(data));
    });
    data.data()(0) = complex(2e-12,1e-12);
    data.data2()(0) = aa::complex_op<double>(4e-24,1e-24,1e-24,1e-24);
    auto scaled = aa::ratio_real_imag(aa::var_result<complex,aa::elliptic_var>(data));
    require(scaled.mean()(0) == 2. && std::abs(scaled.var()(0)-4.) < 1e-14,
            "small denominator changed the componentwise ratio propagation scale");
    data.data()(0) = complex(1e154,1e154);
    data.data2()(0) = aa::complex_op<double>(1e308,1e308,1e308,1e308);
    auto cancelling = aa::ratio_real_imag(aa::var_result<complex,aa::elliptic_var>(data));
    require(cancelling.mean()(0) == 1. && cancelling.var()(0) == 0.
         && cancelling.stderror()(0) == 0., "finite correlated covariance overflowed before cancellation");
    data.data()(0) = complex(1e200,1.);
    data.data2()(0) = aa::complex_op<double>(0.,0.,0.,0.);
    auto large = aa::ratio_real_imag(aa::var_result<complex,aa::elliptic_var>(data));
    require(large.mean()(0) == 1e200 && large.var()(0) == 0. && large.stderror()(0) == 0.,
            "zero covariance acquired overflow from a large finite ratio");
    data.data()(0) = complex(std::numeric_limits<double>::max(),.5);
    throws_estimate<std::overflow_error>([&] {
        aa::ratio_real_imag(aa::var_result<complex,aa::elliptic_var>(data));
    });
    data.data()(0) = complex(.5,.25);
    data.data2()(0).rere() = std::numeric_limits<double>::max()/2;
    throws_estimate<std::overflow_error>([&] {
        aa::ratio_real_imag(aa::var_result<complex,aa::elliptic_var>(data));
    });
    data.data()(0) = complex(2.,1.);
    data.data2()(0).rere() = std::numeric_limits<double>::infinity();
    throws_estimate<std::domain_error>([&] {
        aa::ratio_real_imag(aa::var_result<complex,aa::elliptic_var>(data));
    });
    data.data2()(0).rere() = -1.;
    throws_estimate<std::domain_error>([&] {
        aa::ratio_real_imag(aa::var_result<complex,aa::elliptic_var>(data));
    });
    throws_estimate<aa::finalized_accumulator>([&] {
        aa::ratio_real_imag(aa::var_result<complex,aa::elliptic_var>());
    });
    aa::result ellipse(pooled);
    require(ellipse.count2() == pooled.count2(), "elliptic wrapper lost sample weights");
    throws_estimate<aa::estimate_type_mismatch>([&] { ellipse.stderror<complex>(); });
}
template<class R> void wrapped_result(R const& value) {
    using T = typename aa::traits<R>::value_type;
    aa::result wrapped(value);
    require(wrapped.count() == value.count() && wrapped.count2() == value.count2()
         && wrapped.mean<T>() == value.mean() && wrapped.stderror<T>().isApprox(value.stderror(),1e-13),
            "native result wrapper changed its typed estimate");
}
template<class T> void result_accessors() {
    aa::mean_acc<T> mean(2);
    aa::var_acc<T> variance(2,3);
    aa::cov_acc<T> covariance(2,3);
    aa::batch_acc<T> batch(2,8,3);
    aa::autocorr_acc<T> autocorr(2,3);
    for (size_t i=0; i<53; ++i) {
        auto value = reduction_sample<T>(i,0);
        mean << value; variance << value; covariance << value; batch << value; autocorr << value;
    }
    wrapped_result(variance.result()); wrapped_result(covariance.result());
    wrapped_result(batch.result()); wrapped_result(autocorr.result());
    aa::result mean_only(mean.result());
    throws_estimate<aa::estimate_unavailable>([&] { mean_only.count2(); });
    throws_estimate<aa::estimate_unavailable>([&] { mean_only.stderror<T>(); });
    using other_type = std::conditional_t<std::is_same_v<T,double>,std::complex<double>,double>;
    aa::result wrong_type(variance.result());
    throws_estimate<aa::estimate_type_mismatch>([&] { wrong_type.stderror<other_type>(); });
    aa::result empty;
    throws_estimate<aa::finalized_accumulator>([&] { empty.count2(); });
    throws_estimate<aa::finalized_accumulator>([&] { empty.stderror<T>(); });
}
void signed_statistics(std::string const& filename) {
    // The observable and sign share every sample and every batch boundary.
    aa::scalar_binary_transformer<double> ratio([](double numerator, double sign) {
        if (sign == 0.) throw std::domain_error("undefined ratio: zero mean sign");
        return numerator/sign;
    });
    auto sample = [](size_t i) {
        auto sign = i%5 == 0 ? -1. : 1.;
        return aa::column<double>{sign*(1. + .125*(i%7) + .0625*(i%3)), sign};
    };
    aa::batch_acc<double> uninterrupted(2,8,3), stopped(2,8,3), resumed(1,2,1);
    for (size_t i=0; i<173; ++i) {
        uninterrupted << sample(i);
        if (i<71) stopped << sample(i);
    }
    {
        alps::hdf5::archive ar(filename, "a");
        aa::hdf5_serializer codec(ar, "/signed");
        serialize(codec, "checkpoint", stopped);
    }
    {
        alps::hdf5::archive ar(filename, "r");
        aa::hdf5_serializer codec(ar, "/signed");
        deserialize(codec, "checkpoint", resumed);
    }
    for (size_t i=71; i<173; ++i) resumed << sample(i);
    require(uninterrupted.result() == resumed.result(), "signed joint restart changed batches");
    require(uninterrupted.offset() == resumed.offset(), "signed joint restart changed batch ordering");
    auto joint = roundtrip(filename, "signed-joint", resumed.result());
    auto oracle = signed_oracle(joint.store());
    Eigen::Vector2d gradient(1./oracle.mean(1), -oracle.mean(0)/(oracle.mean(1)*oracle.mean(1)));
    auto variance = double((gradient.transpose()*oracle.covariance*gradient)(0,0));
    auto error = std::sqrt(variance*oracle.count2/(oracle.count*oracle.count));
    auto propagated = aa::transform(aa::linear_prop(), ratio, joint);
    require(std::abs(propagated.mean()(0) - oracle.mean(0)/oracle.mean(1)) < 1e-14,
            "signed ratio mean changed");
    require(std::abs(propagated.var()(0) - variance) < 1e-8*variance,
            "signed ratio discarded joint covariance");
    require(std::abs(propagated.stderror()(0) - error) < 1e-8*error,
            "signed ratio uncertainty disagrees with covariance oracle");
    auto diagonal_only = gradient(0)*gradient(0)*oracle.covariance(0,0)
                       + gradient(1)*gradient(1)*oracle.covariance(1,1);
    require(diagonal_only > 2*variance, "fixture does not distinguish independent from joint errors");
    roundtrip(filename, "signed-linear", propagated);

    // Independent weighted delete-one-bin pseudovalues and uncertainty.
    aa::batch_data<double> pseudovalues(1, joint.num_batches());
    aa::column<double> total = joint.store().batch().rowwise().sum();
    for (size_t i=0; i<joint.num_batches(); ++i) {
        auto weight = double(joint.store().count()(i));
        pseudovalues.count()(i) = joint.store().count()(i);
        if (!weight) continue;
        aa::column<double> leaveout = (total-joint.store().batch().col(i))/(oracle.count-weight);
        pseudovalues.batch()(0,i) = oracle.count*(oracle.mean(0)/oracle.mean(1))
            - (oracle.count-weight)*(leaveout(0)/leaveout(1));
    }
    auto jack_oracle = signed_oracle(pseudovalues);
    auto jack = aa::transform(aa::jackknife_prop(), ratio, joint);
    require(std::abs(jack.mean()(0) - jack_oracle.mean(0)) < 1e-12,
            "signed jackknife mean disagrees with weighted oracle");
    require(std::abs(jack.var()(0) - jack_oracle.covariance(0,0)) < 1e-10,
            "signed jackknife variance disagrees with weighted oracle");
    auto jack_error = std::sqrt(jack_oracle.covariance(0,0)*jack_oracle.count2
                               /(jack_oracle.count*jack_oracle.count));
    require(std::abs(jack.stderror()(0) - jack_error) < 1e-10,
            "signed jackknife error disagrees with weighted oracle");
    require(jack.store().count() == joint.store().count(), "jackknife changed bin weights");
    roundtrip(filename, "signed-jackknife", jack);

    // A fluctuating sign does not create uncertainty in an exactly constant ratio.
    aa::batch_acc<double> constant_ratio(2,8,2), deterministic(2,8,1);
    for (size_t i=0; i<193; ++i) {
        auto sign = i%4 == 0 ? -1. : 1.;
        constant_ratio << aa::column<double>{2*sign, sign};
        deterministic << aa::column<double>{3.,1.};
    }
    auto constant_linear = aa::transform(aa::linear_prop(), ratio, constant_ratio.result());
    auto constant_jack = aa::transform(aa::jackknife_prop(), ratio, constant_ratio.result());
    require(std::isfinite(constant_linear.stderror()(0)) && constant_linear.stderror()(0) < 1e-7,
            "correlated constant ratio acquired uncertainty");
    require(constant_jack.mean()(0) == 2. && constant_jack.stderror()(0) == 0.,
            "constant ratio jackknife lost exact cancellation");
    auto exact = aa::transform(aa::linear_prop(), ratio, deterministic.result());
    require(exact.mean()(0) == 3. && exact.stderror()(0) == 0.,
            "zero input uncertainty produced a nonfinite derivative");

    // Slight negative cancellation is corrected; genuinely invalid variance remains.
    aa::cov_result<double> cancellation(aa::cov_data<double>(2));
    cancellation.store().data() = aa::column<double>{.125,.1};
    cancellation.store().data2() << 1.5625,1.25,1.25,1.;
    cancellation.store().data2() *= .125;
    cancellation.store().count() = 64;
    cancellation.store().count2() = 64.;
    auto cancelled = aa::transform(aa::linear_prop(), ratio, cancellation);
    require(std::isfinite(cancelled.stderror()(0)) && cancelled.stderror()(0) < 1e-7,
            "roundoff in cancelled covariance produced a NaN error");
    cancellation.store().data2() = -Eigen::Matrix2d::Identity();
    auto invalid = aa::transform(aa::linear_prop(), ratio, cancellation);
    require(invalid.var()(0) < 0., "propagation silently zeroed genuinely invalid covariance");

    // Nonzero signs have no arbitrary cutoff, even when their scale is tiny.
    auto derivative = aa::jacobian(ratio, aa::column<double>{2.5e-12,1e-12}, 0.);
    require(std::abs(derivative(0,0)/1e12 - 1.) < 1e-8
         && std::abs(derivative(0,1)/-2.5e12 - 1.) < 1e-8,
            "Jacobian step destroyed near-zero denominator accuracy");
    rejects([&] { aa::linear_prop bad(-1.); });
    rejects([&] { aa::linear_prop bad(std::numeric_limits<double>::quiet_NaN()); });
    rejects([&] { aa::jacobian(ratio, aa::column<double>{1.}, 0.); });
    aa::batch_acc<double> zero_sign(2,8,1), singular_leaveout(2,8,1), insufficient(2,8,3);
    zero_sign << aa::column<double>{1.,1.} << aa::column<double>{1.,-1.};
    rejects([&] { aa::transform(aa::linear_prop(), ratio, zero_sign.result()); });
    rejects([&] { aa::transform(aa::jackknife_prop(), ratio, zero_sign.result()); });
    singular_leaveout << aa::column<double>{1.,1.} << aa::column<double>{1.,-1.}
                     << aa::column<double>{1.,1.};
    rejects([&] { aa::transform(aa::jackknife_prop(), ratio, singular_leaveout.result()); });
    rejects([&] { aa::transform(aa::jackknife_prop(), ratio, insufficient.result()); });
    insufficient << aa::column<double>{1.,1.} << aa::column<double>{2.,1.};
    rejects([&] { aa::transform(aa::jackknife_prop(), ratio, insufficient.result()); });

    // UINT64 counts above PTRDIFF_MAX retain their positive interpretation.
    aa::batch_data<double> huge(1,2);
    huge.count()(0) = uint64_t(1)<<63;
    huge.count()(1) = (uint64_t(1)<<63)-1;
    for (size_t i=0; i<2; ++i) huge.batch()(0,i) = double(huge.count()(i));
    aa::linear_transformer<double> identity(Eigen::Matrix<double,1,1>::Identity());
    auto huge_pseudo = aa::jackknife(huge, identity);
    require(huge_pseudo.count() == huge.count() && huge_pseudo.batch() == huge.batch(),
            "jackknife interpreted uint64 count as negative or lost its small-bin term");
    huge.count()(0) = (uint64_t(1)<<63)-2;
    huge.count()(1) = 1;
    huge.batch()(0,0) = double(huge.count()(0));
    huge.batch()(0,1) = 2.;
    huge_pseudo = aa::jackknife(huge, identity);
    require(huge_pseudo.batch() == huge.batch(), "linear jackknife lost unequal small-bin data");
    aa::batch_data<double> malformed(1,3);
    malformed.count().resize(2);
    malformed.count() << 1,1;
    rejects([&] { aa::jackknife(malformed, identity); });
    aa::batch_data<double> invalid_empty(1,3);
    invalid_empty.count() << 1,1,0;
    invalid_empty.batch() << 1.,2.,3.;
    rejects([&] { aa::jackknife(invalid_empty, identity); });
    huge.count()(0) = std::numeric_limits<uint64_t>::max();
    huge.count()(1) = 1;
    rejects([&] { aa::jackknife(huge, identity); });
}
void eigen_orientation(std::string const& filename) {
    Eigen::Matrix<double,2,3,Eigen::RowMajor> original;
    original << 1., 2., 3., 4., 5., 6.;
    {
        alps::hdf5::archive ar(filename, "a");
        aa::hdf5_serializer codec(ar, "/eigen");
        alps::serialization::serialize(codec, "matrix", original);
        require(ar.extent("/eigen/matrix") == std::vector<size_t>({3,2}), "Eigen physical axes changed");
        std::vector<double> flat(6);
        ar.read("/eigen/matrix", flat.data(), {3,2});
        require(flat == std::vector<double>({1.,4.,2.,5.,3.,6.}), "row-major fallthrough corrupted column-major payload");
        alps::serialization::serialize(codec, "complex-scalar", std::complex<double>(2., -3.));
        require(ar.extent("/eigen/complex-scalar").empty(), "complex scalar acquired an extra axis");
    }
    Eigen::Matrix<double,2,3,Eigen::RowMajor> restored;
    {
        alps::hdf5::archive ar(filename, "r");
        aa::hdf5_serializer codec(ar, "/eigen");
        alps::serialization::deserialize(codec, "matrix", restored);
        std::complex<double> scalar;
        alps::serialization::deserialize(codec, "complex-scalar", scalar);
        require(scalar == std::complex<double>(2., -3.), "complex scalar did not round trip");
    }
    require(original == restored, "row-major matrix did not round trip");
}
void failed_loads(std::string const& filename) {
    aa::cov_acc<double> source(2);
    for (int i=0; i<16; ++i) source << std::vector<double>{double(i), double(2*i)};
    auto original = source.result();
    roundtrip(filename, "bad-covariance", original);
    {
        alps::hdf5::archive ar(filename, "a");
        ar.delete_data("/results/bad-covariance/cov");
        ar.write("/results/bad-covariance/cov", 1.);
    }
    auto target = original;
    {
        alps::hdf5::archive ar(filename, "r");
        aa::hdf5_serializer codec(ar, "/results");
        rejects([&] { deserialize(codec, "bad-covariance", target); });
        require(target == original, "failed result load mutated the destination");
        deserialize(codec, "correlated", target); // traversal recovered after nested failure
    }
    aa::batch_acc<double> acc(2,8,3);
    for (size_t i=0; i<137; ++i) acc << sample<double>(i,2);
    auto before = acc.result();
    auto before_offsets = acc.offset().eval();
    auto before_cursor = acc.cursor();
    {
        alps::hdf5::archive ar(filename, "a");
        aa::hdf5_serializer codec(ar, "/");
        serialize(codec, "bad-state", acc);
        ar.write("/bad-state/cursor/level_position", uint64_t(999));
    }
    {
        alps::hdf5::archive ar(filename, "r");
        aa::hdf5_serializer codec(ar, "/");
        rejects([&] { deserialize(codec, "bad-state", acc); });
    }
    require(before == acc.result(), "failed accumulator load mutated the destination");
    require(before.store().count() == acc.result().store().count()
         && before_offsets == acc.offset()
         && before_cursor.current() == acc.cursor().current()
         && before_cursor.cycle() == acc.cursor().cycle()
         && before_cursor.level() == acc.cursor().level(),
            "failed accumulator load changed continuation state");
    aa::batch_acc<double> sparse(2,8,3);
    sparse << sample<double>(0,2);
    for (int corruption : {0,1,2,3}) {
        auto const& source = corruption == 3 ? sparse : acc;
        {
            alps::hdf5::archive ar(filename, "a");
            aa::hdf5_serializer codec(ar, "/");
            serialize(codec, "bad-layout", source);
            if (corruption == 0) {
                auto offset = source.offset().eval();
                ++offset(0);
                ar.write("/bad-layout/batch/offset", offset.data(), {8});
            } else if (corruption == 1) {
                auto count = source.store().count().eval();
                size_t different = 1;
                while (count(different) == count(0)) ++different;
                std::swap(count(0), count(different)); // Same total, wrong native layout.
                ar.write("/bad-layout/batch/count", count.data(), {8});
            } else if (corruption == 2) {
                ar.write("/bad-layout/cursor/level_position", uint64_t(2)); // In bounds, inconsistent with counts.
            } else {
                auto sums = source.store().batch().eval();
                sums(0,1) = 1.; // Unoccupied slot must remain zero.
                ar.write("/bad-layout/batch/sum", sums.data(), {8,2});
            }
        }
        {
            alps::hdf5::archive ar(filename, "r");
            aa::hdf5_serializer codec(ar, "/");
            rejects([&] { deserialize(codec, "bad-layout", acc); });
        }
        require(before == acc.result() && before_offsets == acc.offset()
             && before_cursor.current() == acc.cursor().current()
             && before_cursor.cycle() == acc.cursor().cycle()
             && before_cursor.level() == acc.cursor().level(),
                "corrupt checkpoint layout changed continuation state");
    }
    aa::batch_acc<double> initial_merge(2,8,3);
    for (size_t i=0; i<25; ++i) initial_merge << sample<double>(i,2);
    {
        alps::hdf5::archive ar(filename, "a");
        aa::hdf5_serializer codec(ar, "/");
        serialize(codec, "inconsistent-state", initial_merge);
        ar.write("/inconsistent-state/cursor/level", uint64_t(64));
    }
    {
        alps::hdf5::archive ar(filename, "r");
        aa::hdf5_serializer codec(ar, "/");
        rejects([&] { deserialize(codec, "inconsistent-state", acc); });
        require(before == acc.result(), "inconsistent cursor published accumulator state");
    }
    {
        alps::hdf5::archive ar(filename, "a");
        ar.write("/results/variance/count2", std::numeric_limits<double>::quiet_NaN());
        ar.write("/results/covariance/count2", 0.);
        ar.write("/results/mean/@size", uint64_t(0));
    }
    {
        alps::hdf5::archive ar(filename, "r");
        aa::hdf5_serializer codec(ar, "/results");
        aa::var_result<double> variance;
        aa::mean_result<double> mean;
        auto covariance = target;
        rejects([&] { deserialize(codec, "variance", variance); });
        rejects([&] { deserialize(codec, "mean", mean); });
        rejects([&] { deserialize(codec, "covariance", covariance); });
        require(!variance.valid() && !mean.valid(), "invalid dimensions/counts published staged results");
        require(covariance == target, "inconsistent squared-weight count replaced existing covariance");
    }
    aa::batch_acc<double> derived(2,8,1);
    for (size_t i=0; i<64; ++i) derived << sample<double>(i,2);
    auto result = roundtrip(filename, "bad-derived", derived.result());
    auto previous = result;
    {
        alps::hdf5::archive ar(filename, "a");
        std::vector<int32_t> wrong_type(2, 0);
        ar.write("/results/bad-derived/mean/error", wrong_type.data(), {2});
    }
    {
        alps::hdf5::archive ar(filename, "r");
        aa::hdf5_serializer codec(ar, "/results");
        rejects([&] { deserialize(codec, "bad-derived", result); });
        require(result == previous, "derived-field failure mutated the destination");
    }
}
}
int main() {
    auto const filename = boost::filesystem::unique_path("statistics.%%%%-%%%%.h5").string();
    try {
        for (auto split : {size_t(0), size_t(1), size_t(3), size_t(8), size_t(23), size_t(24),
                           size_t(25), size_t(64), size_t(137)}) {
            for (uint64_t base : {uint64_t(1), uint64_t(3)}) {
                resume<double>(filename,1,split,base);
                resume<double>(filename,3,split,base);
                resume<std::complex<double>>(filename,2,split,base);
            }
        }
        results(filename);
        mean_tests();
        unequal_merge(filename);
        independent_reductions<double>();
        independent_reductions<std::complex<double>>();
        empty_autocorrelation_reduction<double>();
        empty_autocorrelation_reduction<std::complex<double>>();
        reduction_limits();
        elliptic_reduction();
        batch_reset_and_equality();
        large_batch_weights();
        wrong_sized_append();
        covariance_transform(filename);
        elliptic_signed_ratios(filename);
        result_accessors<double>();
        result_accessors<std::complex<double>>();
        signed_statistics(filename);
        eigen_orientation(filename);
        failed_loads(filename);
        boost::filesystem::remove(filename);
        return 0;
    } catch (std::exception const& error) {
        std::cerr << error.what() << '\n';
        boost::filesystem::remove(filename);
        return 1;
    }
}
