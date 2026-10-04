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
template<class T> void resume(std::string const& filename, size_t components, size_t boundary) {
    aa::batch_acc<T> uninterrupted(components, 8, 3);
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
        for (auto split : {size_t(0), size_t(23), size_t(137)}) {
            resume<double>(filename,1,split);
            resume<double>(filename,3,split);
            resume<std::complex<double>>(filename,2,split);
        }
        results(filename);
        unequal_merge(filename);
        independent_reductions<double>();
        independent_reductions<std::complex<double>>();
        reduction_limits();
        elliptic_reduction();
        batch_reset_and_equality();
        large_batch_weights();
        wrong_sized_append();
        covariance_transform(filename);
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
