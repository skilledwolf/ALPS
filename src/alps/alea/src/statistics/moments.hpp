// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#pragma once
#include <alps/alea/internal/util.hpp>

namespace alps::alea::internal {

// Weighted Chan/Welford update. data is the mean, data2 the centered sum.
// Keeping these centered also makes constant streams exactly noiseless.
template<class Data, class Mean, class Square>
void add_mean(Data& into, Mean const& mean, uint64_t weight, double weight2,
              Square square) {
    if (!weight) return;
    if (weight > std::numeric_limits<uint64_t>::max() - into.count()
        || !std::isfinite(into.count2() + weight2))
        throw std::overflow_error("ALEA moment weight overflow");
    if (!into.count()) into.data() = mean;
    else {
        auto delta = (mean - into.data()).eval();
        double fraction = double(weight) / (into.count() + weight);
        into.data2() += (double(into.count()) * fraction) * square(delta);
        into.data() += fraction * delta;
    }
    into.count() += weight;
    into.count2() += weight2;
}

template<class Data, class Square>
void merge_moments(Data& into, Data const& from, Square square) {
    add_mean(into, from.data(), from.count(), from.count2(), square);
    if (from.count()) into.data2() += from.data2();
}

// Gather only means and weights; sum the already centered moments directly.
// Extra storage is O(ranks * components), even for a covariance matrix. This
// uses the existing sum-only transport and works with deferred reductions.
template<class Data, class Square>
void reduce_moments(Data& data, reducer const& r, Square square) {
    using T = typename Data::value_type;
    using M = typename std::decay_t<decltype(data.data2())>::Scalar;
    auto setup = r.get_setup();
    auto limit = uint64_t(std::numeric_limits<Eigen::Index>::max());
    if (setup.count > limit / data.size()) throw size_mismatch();
    if (r.get_max(!valid_weight_count(data.count(), data.count2())))
        throw std::runtime_error("invalid ALEA squared-weight count");
    typename eigen<T>::matrix means(data.size(), setup.count);
    column<uint64_t> weights(setup.count);
    means.setZero(); weights.setZero();
    if (data.count()) means.col(setup.pos) = data.data();
    weights(setup.pos) = data.count();
    data.unnormalize();
    r.reduce(view<T>(means.data(), means.size()));
    r.reduce(view<uint64_t>(weights.data(), weights.size()));
    r.reduce(view<M>(data.data2().data(), data.data2().size()));
    r.reduce(view<double>(&data.count2(), 1));
    r.commit();
    uint64_t total = 0;
    bool overflow = false;
    // Index the weights: Eigen 3.3 vectors have no begin()/end().
    if (setup.have_result) for (Eigen::Index rank = 0; rank < weights.size(); ++rank) {
        overflow = overflow || weights(rank) > std::numeric_limits<uint64_t>::max() - total;
        total += weights(rank);
    }
    if (r.get_max(setup.have_result && (overflow || !valid_weight_count(total, data.count2()))))
        throw std::overflow_error("ALEA moment weight overflow");
    if (!setup.have_result) return;
    data.count() = 0;
    data.data().setZero();
    for (size_t rank = 0; rank < setup.count; ++rank)
        add_mean(data, means.col(rank), weights(rank), 0., square);
    data.normalize();
}

}
