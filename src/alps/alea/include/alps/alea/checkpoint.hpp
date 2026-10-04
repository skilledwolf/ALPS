// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#pragma once
#include <alps/alea/batch.hpp>
#include <alps/alea/serialize.hpp>
#include <alps/alea/internal/util.hpp>
#include <limits>

namespace alps::alea::internal {
// Only batch_acc currently has a resumable-state codec. Result codecs cannot
// substitute for accumulator checkpoints: the merge cursor and offsets matter.
struct batch_checkpoint {
    // slots+1 is at most 2^63, so each modular addition fits uint64_t.
    static uint64_t multiply_mod(uint64_t a, uint64_t b, uint64_t modulus) {
        uint64_t product = 0;
        a %= modulus;
        for (; b; b >>= 1) {
            if (b & 1) product = (product+a)%modulus;
            a = (a+a)%modulus;
        }
        return product;
    }
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
        if (level > 63 || position >= (level ? shape[0]/2 : shape[0]))
            throw std::runtime_error("invalid ALEA batch checkpoint cursor");
        staged.cursor_.level_ = level;
        staged.cursor_.level_pos_ = position;
        staged.cursor_.factor_ = uint64_t(1)<<level;
        staged.cursor_.skip_ = level ? staged.cursor_.factor_/2 : 0;
        staged.cursor_.current_ = position;
        staged.cursor_.cycle_ = 0;
        if (level) {
            auto const modulus = uint64_t(shape[0])+1;
            // At a level's start current=skip-1; each position advances 2*skip.
            auto const product = multiply_mod(staged.cursor_.skip_, 1+2*position, modulus);
            staged.cursor_.current_ = (product+modulus-1)%modulus;
            uint64_t power = 1;
            for (uint64_t i=0; i<level; ++i) {
                if (power == 1) ++staged.cursor_.cycle_;
                power = (power+power)%modulus;
            }
        }
        if (base_size > std::numeric_limits<uint64_t>::max()/staged.cursor_.factor_)
            throw std::overflow_error("ALEA current batch size overflows");
        auto target = base_size * staged.cursor_.factor_;
        uint64_t count = 0;
        for (auto n : staged.store_->count()) {
            if (n > target || count > std::numeric_limits<uint64_t>::max()-n)
                throw std::runtime_error("invalid ALEA batch checkpoint count");
            count += n;
        }
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
