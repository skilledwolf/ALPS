// Copyright (C) 2026 ALPS Collaboration.
// SPDX-License-Identifier: MIT

#include <alps/alea.h>
#include <alps/hdf5.hpp>

#include <filesystem>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <type_traits>

namespace {
void require(bool condition, char const* message) {
    if (!condition) throw std::runtime_error(message);
}

template<typename T> bool equal(T const& lhs, T const& rhs) {
    if constexpr (std::is_arithmetic_v<T>) return lhs == rhs;
    else {
        if (lhs.size() != rhs.size()) return false;
        for (std::size_t i = 0; i < lhs.size(); ++i)
            if (!equal(lhs[i], rhs[i])) return false;
        return true;
    }
}

template<typename T> T sample(std::size_t i) {
    double value = double(i % 17) - 8.0;
    if constexpr (std::is_arithmetic_v<T>) return value;
    else return T{value, 2.0 * value + 3.0};
}

template<typename Observable> void observable_roundtrip(alps::hdf5::archive& ar,
                                                        std::string const& path) {
    using T = typename Observable::value_type;
    using Data = alps::SimpleObservableData<T>;
    Observable observable("energy");
    for (std::size_t i = 0; i < 256; ++i) observable << sample<T>(i);
    Data original(observable);
    ar << alps::make_pvp(path, original);
    Data restored;
    ar >> alps::make_pvp(path, restored);
    require(restored.count() == original.count(), "observable count was not restored");
    require(restored.bin_size() == original.bin_size(), "observable bin size was not restored");
    require(equal(restored.mean(), original.mean()), "observable mean changed");
    require(equal(restored.error(), original.error()), "observable error changed");
    require(equal(restored.variance(), original.variance()), "observable variance changed");
    require(equal(restored.tau(), original.tau()), "observable tau changed");
    require(equal(restored.bins(), original.bins()), "observable bin sums changed");
    require(equal(restored.any_converged_errors(), original.any_converged_errors()),
            "observable convergence metadata changed");
    uint64_t maximum = 0;
    ar >> alps::make_pvp(path + "/timeseries/data2/@maxbinnum", maximum);
    require(maximum == observable.max_bin_number(), "squared-bin metadata was not written");

    // Released evaluators omitted binsize; infer it from count/bin count.
    ar.delete_attribute(path + "/timeseries/data/@binsize");
    ar >> alps::make_pvp(path, restored);
    require(restored.bin_size() == original.bin_size(), "official evaluator bin-size inference failed");

    auto before = restored;
    ar << alps::make_pvp(path + "/count", uint64_t(17));
    ar.delete_data(path + "/jacknife/data");
    ar << alps::make_pvp(path + "/jacknife/data", std::string("malformed numeric state"));
    bool rejected = false;
    try { ar >> alps::make_pvp(path, restored); }
    catch (std::exception const&) { rejected = true; }
    require(rejected, "malformed late observable field was accepted");
    require(restored.count() == before.count() && restored.bin_size() == before.bin_size(),
            "failed observable load changed count or bin size");
    require(equal(restored.mean(), before.mean()) && equal(restored.error(), before.error()),
            "failed observable load changed statistics");
    require(equal(restored.bins(), before.bins()) && equal(restored.tau(), before.tau()),
            "failed observable load changed optional state");
    ar << alps::make_pvp(path, original);
    ar.delete_data(path + "/tau/value");
    ar << alps::make_pvp(path + "/tau/value/0", std::string("not numeric"));
    rejected = false;
    try { ar >> alps::make_pvp(path, restored); }
    catch (std::exception const&) { rejected = true; }
    require(rejected && equal(restored.tau(), before.tau()),
            "malformed optional group was ignored or changed state");
    ar.delete_group(path + "/tau/value");
    ar << alps::make_pvp(path, original);

    // Absent optionals replace existing values, and saving removes only owned leaves.
    ar << alps::make_pvp(path + "/variance/application", 19);
    ar << alps::make_pvp(path + "/tau/application", 23);
    ar << alps::make_pvp(path + "/jacknife/application", 29);
    ar.delete_data(path + "/variance/value");
    ar.delete_data(path + "/tau/value");
    ar.delete_data(path + "/jacknife/data");
    ar >> alps::make_pvp(path, restored);
    require(!restored.has_variance() && !restored.has_tau(), "absent observable statistics survived load");
    ar << alps::make_pvp(path, original);
    ar << alps::make_pvp(path, restored);
    for (auto field : {"/variance/value", "/tau/value", "/jacknife/data"})
        require(!ar.is_data(path + field), "obsolete observable optional leaf survived save");
    for (auto field : {"/variance/application", "/tau/application", "/jacknife/application"})
        require(ar.is_data(path + field), "observable save removed unrelated application data");

    ar << alps::make_pvp(path, Data());
    ar >> alps::make_pvp(path, restored);
    require(restored.count() == 0 && restored.bin_number() == 0 && restored.bin_size() == 0,
            "empty observable did not replace populated state");
    require(!restored.has_variance() && !restored.has_tau(), "empty observable retained statistics");
    for (auto field : {"/variance/value", "/tau/value"})
        require(!ar.is_data(path + field), "empty observable save retained optional leaf");

    // A result discriminator left in a reused group must not disable sum decoding.
    ar << alps::make_pvp(path + "/@cannotrebin", true);
    ar << alps::make_pvp(path, original);
    require(!ar.is_attribute(path + "/@cannotrebin"), "observable retained result discriminator");
    // Rebinning a loaded evaluator needs its restored bin size.
    ar >> alps::make_pvp(path, restored);
    auto rebinned = original;
    restored.set_bin_size(original.bin_size() * 2);
    rebinned.set_bin_size(original.bin_size() * 2);
    require(equal(restored.mean(), rebinned.mean()) && equal(restored.error(), rebinned.error()),
            "loaded observable rebinning changed statistics");

    Observable single_observable("energy");
    single_observable << sample<T>(0);
    ar << alps::make_pvp(path + "/single_sample", single_observable);
    require(!ar.is_data(path + "/single_sample/mean/error"),
            "raw first-sample observable unexpectedly serialized an error");
}

void assignment_preserves_evaluation_method() {
    std::istringstream xml("<SCALAR_AVERAGE><COUNT>1</COUNT><MEAN>3</MEAN>"
                           "<ERROR method=\"bootstrap\" converged=\"yes\">0.25</ERROR>"
                           "</SCALAR_AVERAGE>");
    auto tag = alps::parse_tag(xml);
    std::string label;
    alps::SimpleObservableData<double> source(xml, tag, label), assigned;
    assigned = source;
    require(assigned.evaluation_method(alps::Error) == "bootstrap",
            "observable assignment dropped evaluation method");
    assigned = assigned;
    require(assigned.mean() == 3.0 && assigned.evaluation_method(alps::Error) == "bootstrap",
            "observable self-assignment changed state");
}
} // namespace

int main() {
    const std::string filename = "alea-checkpoint-contracts.h5";
    try {
        assignment_preserves_evaluation_method();
        {
            alps::hdf5::archive ar(filename, "w");
            observable_roundtrip<alps::RealObservable>(ar, "/scalar");
            observable_roundtrip<alps::RealVectorObservable>(ar, "/vector");
        }
        std::filesystem::remove(filename);
    } catch (std::exception const& error) {
        std::filesystem::remove(filename);
        std::cerr << error.what() << '\n';
        return 1;
    }
}
