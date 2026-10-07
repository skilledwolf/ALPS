// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#pragma once
#include <alps/alea/batch.hpp>
#include <alps/alea/mean.hpp>
#include <alps/alea/covariance.hpp>
#include <alps/alea/autocorr.hpp>
#include <alps/alea/serialize.hpp>
#include <alps/alea/internal/util.hpp>
#include <algorithm>
#include <limits>

namespace alps::alea::internal {
// Result codecs cannot substitute for accumulator checkpoints: partial bins,
// centered sums and the merge cursor must survive without normalization.
struct batch_checkpoint {
    template<class T> static void save(serializer& s, std::string const& key, batch_acc<T> const& acc) {
        check_valid(acc);
        result_sentry group(s, key, result_kind::batch_accumulator);
        using alps::serialization::serialize;
        serialize(s, "base_size", uint64_t(acc.base_size_));
        {
            serializer_sentry batch(s, "batch");
            serialize(s, "sum", acc.store_->batch());
            serialize(s, "count", acc.store_->count());
            serialize(s, "offset", acc.offset_);
        }
        {
            serializer_sentry cursor(s, "cursor");
            serialize(s, "level", uint64_t(acc.cursor_.level_));
            serialize(s, "level_position", uint64_t(acc.cursor_.level_pos_));
        }
    }
    template<class T> static void load(deserializer& s, std::string const& key, batch_acc<T>& acc) {
        result_reader_sentry group(s, key, result_kind::batch_accumulator);
        using alps::serialization::deserialize;
        uint64_t base_size;
        deserialize(s, "base_size", base_size);
        std::vector<size_t> shape;
        { deserializer_sentry batch(s, "batch"); shape = s.get_shape("sum"); }
        if (shape.size() != 2 || !shape[1]
                || shape[0] > size_t(std::numeric_limits<Eigen::Index>::max())
                || shape[1] > size_t(std::numeric_limits<Eigen::Index>::max())) throw size_mismatch();
        if (base_size > std::numeric_limits<size_t>::max()) throw std::overflow_error("ALEA batch size overflows");
        batch_acc<T> staged(shape[1], shape[0], base_size);
        auto expected_offset = staged.offset_.eval();
        typename eigen<uint64_t>::row expected_count = eigen<uint64_t>::row::Zero(shape[0]);
        {
            deserializer_sentry batch(s, "batch");
            deserialize(s, "sum", staged.store_->batch());
            deserialize(s, "count", staged.store_->count());
            deserialize(s, "offset", staged.offset_);
        }
        uint64_t level, position;
        {
            deserializer_sentry cursor(s, "cursor");
            deserialize(s, "level", level);
            deserialize(s, "level_position", position);
        }
        uint64_t count = 0;
        for (size_t i=0; i<shape[0]; ++i) {
            auto n = staged.store_->count()(i);
            if ((!n && !staged.store_->batch().col(i).isZero(0))
                    || count > std::numeric_limits<uint64_t>::max()-n)
                throw std::runtime_error("invalid ALEA batch checkpoint count");
            count += n;
        }
        // Unit samples determine the entire layout. Replay complete batches,
        // including native merges, in O(slots * log(count)) rather than samples.
        uint64_t filled = 0;
        while (filled < count) {
            if (base_size > std::numeric_limits<uint64_t>::max()/staged.cursor_.factor())
                throw std::overflow_error("ALEA current batch size overflows");
            auto added = std::min(base_size * staged.cursor_.factor(), count-filled);
            expected_count(staged.cursor_.current()) += added;
            filled += added;
            if (filled == count) break;
            ++staged.cursor_;
            if (staged.cursor_.merge_mode()) {
                auto into = staged.cursor_.merge_into(), current = staged.cursor_.current();
                expected_count(into) += expected_count(current);
                expected_count(current) = 0;
                expected_offset(into) = std::min(expected_offset(into), expected_offset(current));
                expected_offset(current) = filled;
            }
        }
        if (level != staged.cursor_.level_ || position != staged.cursor_.level_pos_
                || expected_count != staged.store_->count() || expected_offset != staged.offset_)
            throw std::runtime_error("inconsistent ALEA batch checkpoint layout");
        acc = std::move(staged);
    }
};

struct accumulator_checkpoint {
    template<class T> static result_kind kind(mean_acc<T> const&) { return result_kind::mean_accumulator; }
    template<class T, class S> static result_kind kind(var_acc<T,S> const&) { return result_kind::variance_accumulator; }
    template<class T, class S> static result_kind kind(cov_acc<T,S> const&) { return result_kind::covariance_accumulator; }

    template<class A> static constexpr bool mean_only = std::is_same_v<A,mean_acc<typename A::value_type>>;

    template<class A> static void save(serializer& s, std::string const& key, A const& acc) {
        check_valid(acc);
        result_sentry group(s,key,kind(acc));
        using alps::serialization::serialize;
        serialize(s,"value",acc.store_->data());
        serialize(s,"count",acc.store_->count());
        if constexpr (!mean_only<A>) {
            serialize(s,"centered_moment",acc.store_->data2());
            serialize(s,"count2",acc.store_->count2());
            serialize(s,"batch_size",acc.current_.target());
            serialize(s,"partial_sum",acc.current_.sum());
            serialize(s,"partial_count",acc.current_.count());
        }
    }
    template<class A> static void load(deserializer& s, std::string const& key, A& acc) {
        result_reader_sentry group(s,key,kind(acc));
        using alps::serialization::deserialize;
        auto shape=s.get_shape("value");
        if (shape.size()!=1 || !shape[0] || shape[0]>size_t(std::numeric_limits<Eigen::Index>::max()))
            throw size_mismatch();
        A staged(shape[0]);
        auto& data=*staged.store_;
        deserialize(s,"value",data.data());
        deserialize(s,"count",data.count());
        if (!data.count() && !data.data().isZero(0))
            throw std::runtime_error("nonzero empty ALEA accumulator");
        if constexpr (!mean_only<A>) {
            // Compare whole matrices: complex_op moments have no abs() for
            // isZero, and Eigen 3.3 array-scalar == breaks C++20 rewriting.
            using moment = std::decay_t<decltype(data.data2())>;
            deserialize(s,"centered_moment",data.data2());
            deserialize(s,"count2",data.count2());
            deserialize(s,"batch_size",staged.current_.target());
            deserialize(s,"partial_sum",staged.current_.sum());
            deserialize(s,"partial_count",staged.current_.count());
            auto const& partial=staged.current_;
            if (!partial.target() || partial.count()>=partial.target()
                || partial.count()>std::numeric_limits<uint64_t>::max()-data.count()
                || (!partial.count() && !partial.sum().isZero(0))
                || !valid_weight_count(data.count(),data.count2())
                || data.count2()>double(data.count())*data.count()
                || (!data.count() && data.data2()!=moment::Zero(data.data2().rows(),data.data2().cols())))
                throw std::runtime_error("invalid ALEA centered-moment checkpoint");
        }
        acc=std::move(staged);
    }
    template<class T> static void save(serializer& s, std::string const& key, autocorr_acc<T> const& acc) {
        check_valid(acc);
        result_sentry group(s,key,result_kind::autocorr_accumulator);
        using alps::serialization::serialize;
        serialize(s,"size",uint64_t(acc.size_));
        serialize(s,"count",uint64_t(acc.count_));
        serialize(s,"batch_size",uint64_t(acc.batch_size_));
        serialize(s,"granularity",uint64_t(acc.granularity_));
        serializer_sentry levels(s,"levels");
        for (size_t i=0;i<acc.nlevel();++i) save(s,std::to_string(i),acc.level_[i]);
    }
    template<class T> static void load(deserializer& s, std::string const& key, autocorr_acc<T>& acc) {
        result_reader_sentry group(s,key,result_kind::autocorr_accumulator);
        using alps::serialization::deserialize;
        uint64_t size,count,base,granularity;
        deserialize(s,"size",size); deserialize(s,"count",count);
        deserialize(s,"batch_size",base); deserialize(s,"granularity",granularity);
        if (!size || size>uint64_t(std::numeric_limits<Eigen::Index>::max())
            || count>std::numeric_limits<size_t>::max() || base>std::numeric_limits<size_t>::max()
            || granularity>std::numeric_limits<size_t>::max()) throw size_mismatch();
        autocorr_acc<T> staged(size,base,granularity);
        staged.count_=count;
        while (count>=staged.nextlevel_) staged.add_level();
        deserializer_sentry levels(s,"levels");
        for (size_t i=0;i<staged.nlevel();++i) {
            auto width=staged.level_[i].batch_size();
            auto available=i ? count/staged.level_[i-1].batch_size()*staged.level_[i-1].batch_size() : count;
            load(s,std::to_string(i),staged.level_[i]);
            auto const& level=staged.level_[i];
            auto full=available/width*width;
            if (level.size()!=size || level.batch_size()!=width || level.store().count()!=full
                || level.current().count()!=available-full || level.store().count2()!=double(full)*width)
                throw std::runtime_error("inconsistent ALEA autocorrelation checkpoint hierarchy");
        }
        acc=std::move(staged);
    }
};

}
namespace alps::alea {
template<class T> void serialize(serializer& s, std::string const& key, batch_acc<T> const& acc) {
    internal::batch_checkpoint::save(s, key, acc);
}
template<class T> void deserialize(deserializer& s, std::string const& key, batch_acc<T>& acc) {
    internal::batch_checkpoint::load(s, key, acc);
}

template<class A, std::enable_if_t<is_alea_acc<A>::value,int> = 0>
void serialize(serializer& s, std::string const& key, A const& acc) {
    internal::accumulator_checkpoint::save(s,key,acc);
}
template<class A, std::enable_if_t<is_alea_acc<A>::value,int> = 0>
void deserialize(deserializer& s, std::string const& key, A& acc) {
    internal::accumulator_checkpoint::load(s,key,acc);
}
}
