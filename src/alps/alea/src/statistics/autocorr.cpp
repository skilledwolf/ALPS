/*
 * Copyright (C) 1998-2018 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 * For use in publications, see ACKNOWLEDGE.TXT
 */
#include <alps/alea/autocorr.hpp>
#include <limits>
#include <alps/alea/serialize.hpp>

#include <alps/alea/internal/util.hpp>
#include <alps/alea/internal/format.hpp>

namespace alps { namespace alea {

template <typename T>
autocorr_acc<T>::autocorr_acc(size_t size, uint64_t batch_size, size_t granularity)
    : size_(size)
    , batch_size_(batch_size)
    , count_(0)
    , nextlevel_(batch_size)
    , granularity_(granularity)
    , level_()
{
    if (granularity < 2) throw std::invalid_argument("ALEA autocorrelation granularity must be at least two");
    level_.push_back(var_acc<T>(size, batch_size));
}

template <typename T>
void autocorr_acc<T>::reset()
{
    count_ = 0;
    nextlevel_ = batch_size_;
    level_.clear();
    level_.push_back(var_acc<T>(size_, batch_size_));
}

template <typename T>
void autocorr_acc<T>::set_size(size_t size)
{
    size_ = size;
    reset();
}

template <typename T>
void autocorr_acc<T>::set_batch_size(uint64_t batch_size)
{
    // TODO: handle the case where we just discard levels more gracefully
    if (!batch_size) throw std::invalid_argument("ALEA batch size must be positive");
    batch_size_ = batch_size;
    reset();
}

template <typename T>
void autocorr_acc<T>::set_granularity(size_t granularity)
{
    // TODO: handle the case where we just discard levels more gracefully
    if (granularity < 2) throw std::invalid_argument("ALEA autocorrelation granularity must be at least two");
    granularity_ = granularity;
    reset();
}

template <typename T>
void autocorr_acc<T>::add_level()
{
    // add a new level on top and push back the nextlevel
    if (nextlevel_ > std::numeric_limits<size_t>::max()/granularity_)
        throw std::overflow_error("ALEA autocorrelation hierarchy size overflows");
    auto next = nextlevel_ * granularity_;
    level_.push_back(var_acc<T>(size(), next));
    nextlevel_ = next;
}

template <typename T>
void autocorr_acc<T>::add(const computed<T> &source, uint64_t count)
{
    assert(count_ < nextlevel_);
    internal::check_valid(*this);
    if (source.size() != size()) throw size_mismatch();

    // if we require next level, then do it!
    if (count > std::numeric_limits<size_t>::max()-count_)
        throw std::overflow_error("ALEA sample count overflows");
    if (count_ + count >= nextlevel_) add_level();
    count_ += count;

    // now add current element at the bottom and watch it propagate
    level_[0].add(source, count, level_.data() + 1);
}

template <typename T>
autocorr_result<T> autocorr_acc<T>::result() const
{
    internal::check_valid(*this);
    autocorr_result<T> result;
    autocorr_acc<T>(*this).finalize_to(result);
    return result;
}

template <typename T>
autocorr_result<T> autocorr_acc<T>::finalize()
{
    autocorr_result<T> result;
    finalize_to(result);
    return result;
}

template <typename T>
void autocorr_acc<T>::finalize_to(autocorr_result<T> &result)
{
    internal::check_valid(*this);
    result.level_.resize(level_.size());

    // Finalize each level.
    // NOTE: it is imperative to do this in bottom-up order, since it collects
    //       the left-over data in the current batch from the low-lying levels
    //       and propagates them upwards.
    for (size_t i = 0; i != level_.size() - 1; ++i)
        level_[i].finalize_to(result.level_[i], level_.data() + i + 1);

    level_[nlevel() - 1].finalize_to(result.level_[nlevel() - 1], nullptr);
    level_.clear();     // signal invalidity
}

template class autocorr_acc<double>;
template class autocorr_acc<std::complex<double> >;

template <typename T>
bool operator==(const autocorr_result<T> &r1, const autocorr_result<T> &r2)
{
    if(r1.nlevel() != r2.nlevel()) return false;
    for(size_t i = 0; i < r1.nlevel(); ++i) {
        if(r1.level(i) != r2.level(i))
            return false;
    }
    return true;
}

template bool operator==(const autocorr_result<double> &r1,
                         const autocorr_result<double> &r2);
template bool operator==(const autocorr_result<std::complex<double>> &r1,
                         const autocorr_result<std::complex<double>> &r2);

template <typename T>
uint64_t autocorr_result<T>::batch_size(size_t i) const
{
    return level_[i].batch_size();
}

template <typename T>
size_t autocorr_result<T>::find_level(uint64_t min_samples) const
{
    // TODO: this can be done in O(1)
    for (unsigned i = nlevel(); i != 0; --i) {
        if (level(i-1).count() / level(i-1).batch_size() >= min_samples)
            return i - 1;
    }
    return 0;
}

template <typename T>
double autocorr_result<T>::count2() const
{
    size_t lvl = find_level(DEFAULT_MIN_SAMPLES);

    return level_[lvl].count2();
}

template <typename T>
column<typename autocorr_result<T>::var_type> autocorr_result<T>::var() const
{
    size_t lvl = find_level(DEFAULT_MIN_SAMPLES);

    // The question is whether to return the ensemble variance or the batch
    // variance here.  We opt for the batch variance, because it makes
    // everything consistent with each other.
    //return level_[lvl].batch_size() * level_[lvl].var();
    return level_[lvl].var();
}

template <typename T>
column<typename autocorr_result<T>::var_type> autocorr_result<T>::stderror() const
{
    size_t lvl = find_level(DEFAULT_MIN_SAMPLES);

    return level_[lvl].stderror();
}


template <typename T>
double autocorr_result<T>::observations() const
{
    size_t lvl = find_level(DEFAULT_MIN_SAMPLES);

    return level_[lvl].observations();
}

template <typename T>
column<typename autocorr_result<T>::var_type> autocorr_result<T>::tau() const
{
    size_t lvl = find_level(DEFAULT_MIN_SAMPLES);
    const column<var_type> &var0 = level_[0].var();
    const column<var_type> &varn = level_[lvl].var();
    const double n = level_[lvl].batch_size();

    return (0.5 * n * varn.array() / var0.array() - 0.5).matrix();
}

template <typename T>
column<int> autocorr_result<T>::converged_errors() const
{
    internal::check_valid(*this);
    column<int> result = column<int>::Ones(size());
    const size_t last = find_level(DEFAULT_MIN_SAMPLES);
    if (last < 3) return result;
    auto reference = level_[last].stderror();
    auto previous = {level_[last-3].stderror(), level_[last-2].stderror(), level_[last-1].stderror()};
    for (size_t component = 0; component < size(); ++component) {
        int status = 0;
        if (!std::isfinite(reference(component))) continue;
        for (auto const& estimate : previous) {
            auto error = estimate(component);
            if (!std::isfinite(error)) { status = 1; break; }
            // Retain the strongest evidence of a rising error: a subsequent
            // plateau must not erase an earlier failed comparison.
            status = std::max(status, error < 0.824 * reference(component) ? 2 :
                                      error < 0.9 * reference(component) ? 1 : 0);
        }
        result(component) = status;
    }
    return result;
}

template <typename T>
void autocorr_result<T>::reduce(const reducer &r)
{
    bool complete = valid();
    for (auto const& level : level_) complete = complete && level.valid();
    auto setup = internal::check_reduction(r, complete, complete ? size() : 0);
    bool malformed = false;
    for (auto const& level : level_)
        malformed = malformed || level.size() != size() || level.count() != count()
                    || !internal::valid_weight_count(level.count(), level.count2());
    if (r.get_max(malformed)) throw size_mismatch();
    if (r.get_max(nlevel() > size_t(std::numeric_limits<int64_t>::max())))
        throw size_mismatch();
    // Empty runs impose no limit on the common depth. If every run is empty,
    // retain one neutral level; otherwise every retained level has all samples.
    auto depth = r.get_max(count() ? -static_cast<int64_t>(nlevel())
                                  : std::numeric_limits<int64_t>::min());
    auto shared_levels = depth == std::numeric_limits<int64_t>::min() ? 1 : -depth;
    if (shared_levels <= 0 || (count() && uint64_t(shared_levels) > nlevel()))
        throw size_mismatch();
    autocorr_result staged(*this);
    staged.level_.resize(shared_levels, level_result_type(var_data<T>(size())));
    for (auto &level : staged.level_) level.reduce_unchecked(r);
    if (!setup.have_result) staged.level_.clear();
    level_.swap(staged.level_);
}

template class autocorr_result<double>;
template class autocorr_result<std::complex<double> >;


template <typename T>
void serialize(serializer &s, const std::string &key, const autocorr_result<T> &self)
{
    internal::check_valid(self);
    internal::result_sentry group(s, key, internal::result_kind::autocorrelation);

    // Write size and nlevel as 64-bit integers for consistency
    serialize(s, "@size", static_cast<uint64_t>(self.size()));
    serialize(s, "@nlevel", static_cast<uint64_t>(self.nlevel()));

    // One dataset per level field, levels first: a group per level would
    // multiply the HDF5 objects of every observable by its level count.
    typename eigen<uint64_t>::row count(self.nlevel());
    typename eigen<double>::row count2(self.nlevel());
    typename eigen<T>::matrix mean(self.size(), self.nlevel());
    typename eigen<typename autocorr_result<T>::var_type>::matrix var(self.size(), self.nlevel());
    for (size_t i = 0; i != self.nlevel(); ++i) {
        auto const& level = self.level_[i].store();
        count(i) = level.count();
        count2(i) = level.count2();
        mean.col(i) = level.data();
        var.col(i) = level.data2();
    }
    {
    internal::serializer_sentry subgroup(s, "level");
    serialize(s, "count", count);
    serialize(s, "count2", count2);
    serialize(s, "mean", mean);
    serialize(s, "var", var);
    }

    {
    internal::serializer_sentry subgroup(s, "mean");
    serialize(s, "value", self.mean());
    serialize(s, "error", self.stderror());
    }
}

template <typename T>
void deserialize(deserializer &s, const std::string &key, autocorr_result<T> &self)
{
    typedef typename autocorr_result<T>::var_type var_type;
    internal::result_reader_sentry group(s, key, internal::result_kind::autocorrelation);
    autocorr_result<T> staged;

    uint64_t new_size;
    deserialize(s, "@size", new_size);
    if (!new_size || new_size > uint64_t(std::numeric_limits<Eigen::Index>::max())) throw size_mismatch();

    // first deserialize the fundamentals and make sure that the target fits
    uint64_t new_nlevel;
    deserialize(s, "@nlevel", new_nlevel);
    if (!new_nlevel || new_nlevel > std::numeric_limits<size_t>::max()) throw size_mismatch();
    if (new_nlevel > uint64_t(std::numeric_limits<Eigen::Index>::max()) / new_size) throw size_mismatch();
    const auto size = Eigen::Index(new_size), nlevel = Eigen::Index(new_nlevel);
    typename eigen<uint64_t>::row count(nlevel);
    typename eigen<double>::row count2(nlevel);
    typename eigen<T>::matrix mean(size, nlevel);
    typename eigen<var_type>::matrix var(size, nlevel);
    {
    internal::deserializer_sentry subgroup(s, "level");
    const std::vector<size_t> levels{size_t(nlevel)}, fields{size_t(nlevel), size_t(size)};
    if (s.get_shape("count") != levels || s.get_shape("count2") != levels
            || s.get_shape("mean") != fields || s.get_shape("var") != fields)
        throw size_mismatch();
    deserialize(s, "count", count);
    deserialize(s, "count2", count2);
    deserialize(s, "mean", mean);
    deserialize(s, "var", var);
    }
    staged.level_.reserve(new_nlevel);
    for (Eigen::Index i = 0; i != nlevel; ++i) {
        if (!internal::valid_weight_count(count(i), count2(i)))
            throw std::runtime_error("invalid ALEA squared-weight count");
        var_data<T, circular_var> level(new_size);
        level.count() = count(i);
        level.count2() = count2(i);
        level.data() = mean.col(i);
        level.data2() = var.col(i);
        staged.level_.emplace_back(level);
    }
    size_t scalar_size = staged.size();
    {
    internal::deserializer_sentry subgroup(s, "mean");
    s.read("value", ndview<T>(nullptr, &scalar_size, 1)); // discard
    s.read("error", ndview<var_type>(nullptr, &scalar_size, 1)); // discard
    }

    self.level_.swap(staged.level_);
}

template void serialize(serializer &, const std::string &key, const autocorr_result<double> &);
template void serialize(serializer &, const std::string &key, const autocorr_result<std::complex<double>> &);

template void deserialize(deserializer &, const std::string &key, autocorr_result<double> &);
template void deserialize(deserializer &, const std::string &key, autocorr_result<std::complex<double> > &);

template <typename T>
std::ostream &operator<<(std::ostream &str, const autocorr_result<T> &self)
{
    internal::check_valid(self);
    internal::format_sentry sentry(str);
    verbosity verb = internal::get_format(str, PRINT_TERSE);

    if (verb == PRINT_VERBOSE)
        str << "<X> = ";
    str << self.mean() << " +- " << self.stderror();

    if (verb == PRINT_VERBOSE) {
        str << "\nLevels:" << PRINT_TERSE;
        for (const var_result<T> &curr : self.level_)
            str << "\n  " << curr;
    }
    return str;
}

template std::ostream &operator<<(std::ostream &, const autocorr_result<double> &);
template std::ostream &operator<<(std::ostream &, const autocorr_result<std::complex<double>> &);

}}
