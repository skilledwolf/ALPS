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

namespace {
namespace aa = alps::alea;
void require(bool condition, char const* message) {
    if (!condition) throw std::runtime_error(message);
}
template<class F> void rejects(F operation) {
    bool failed = false;
    try { operation(); } catch (std::exception const&) { failed = true; }
    require(failed, "malformed statistics checkpoint was accepted");
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
    require(result.mean()(0) == 2. && result.var()(0) == 2.
         && result.cov()(0,0) == 2. && result.stderror()(0) == 1.,
            "large integer batch weights overflowed squared-weight counts");
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
        batch_reset_and_equality();
        large_batch_weights();
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
