// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#pragma once
#include <alps/alea/batch.hpp>
#include <alps/alea/serialize.hpp>
#include <alps/alea/internal/util.hpp>
#include <algorithm>
#include <limits>

namespace alps::alea::internal {
// Only batch_acc currently has a resumable-state codec. Result codecs cannot
// substitute for accumulator checkpoints: the merge cursor and offsets matter.
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
}
namespace alps::alea {
template<class T> void serialize(serializer& s, std::string const& key, batch_acc<T> const& acc) {
    internal::batch_checkpoint::save(s, key, acc);
}
template<class T> void deserialize(deserializer& s, std::string const& key, batch_acc<T>& acc) {
    internal::batch_checkpoint::load(s, key, acc);
}
}
