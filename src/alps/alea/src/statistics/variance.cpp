/*
 * Copyright (C) 1998-2018 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 * For use in publications, see ACKNOWLEDGE.TXT
 */
#include <alps/alea/variance.hpp>
#include "moments.hpp"
#include <limits>
#include <cmath>
#include <alps/alea/util.hpp>
#include <alps/alea/serialize.hpp>

#include <alps/alea/internal/util.hpp>
#include <alps/alea/internal/format.hpp>

namespace alps { namespace alea {

template <typename T, typename Str>
var_data<T,Str>::var_data(size_t size)
    : data_(size)
    , data2_(size)
{
    reset();
}

template <typename T, typename Str>
void var_data<T,Str>::reset()
{
    data_.fill(0);
    data2_.fill(0);
    count_ = 0;
    count2_ = 0;
}

template <typename T, typename Str>
void var_data<T,Str>::normalize()
{
    // This also works for count_ == 0
    if (!count_) data_.fill(std::numeric_limits<double>::quiet_NaN());

    // In case of zero unbiased information, the variance is infinite.
    // However, data2_ is 0 in this case as well, so we need to handle it
    // specially to avoid 0/0 = nan while propagating intrinsic NaN's.
    const double nunbiased = count_ - count2_/count_;
    if (nunbiased == 0) {
        data2_ = data2_.array().isNaN().select(data2_, INFINITY);
    } else {
        // HACK: this is written in out-of-place notation to work around Eigen
        data2_ = data2_ / nunbiased;
    }
}

template <typename T, typename Str>
void var_data<T,Str>::unnormalize()
{
    // "empty" sets must be handled specially here because of NaNs
    if (count_ == 0) {
        reset();
        return;
    }

    // Care must be taken again for zero unbiased info since inf/0 is NaN.
    const double nunbiased = count_ - count2_/count_;
    if (nunbiased == 0)
        data2_ = data2_.array().isNaN().select(data2_, 0);
    else
        data2_ = data2_ * nunbiased;

}

template class var_data<double>;
template class var_data<std::complex<double>, circular_var>;
template class var_data<std::complex<double>, elliptic_var>;


template <typename T, typename Str>
var_acc<T,Str>::var_acc(size_t size, uint64_t batch_size)
    : store_(new var_data<T,Str>(size))
    , current_(size, batch_size)
{ }

// We need an explicit copy constructor, as we need to copy the data
template <typename T, typename Str>
var_acc<T,Str>::var_acc(const var_acc &other)
    : store_(other.store_ ? new var_data<T,Str>(*other.store_) : nullptr)
    , current_(other.current_)
{ }

template <typename T, typename Str>
var_acc<T,Str> &var_acc<T,Str>::operator=(const var_acc &other)
{
    store_.reset(other.store_ ? new var_data<T,Str>(*other.store_) : nullptr);
    current_ = other.current_;
    return *this;
}

template <typename T, typename Str>
void var_acc<T,Str>::reset()
{
    current_.reset();
    if (valid())
        store_->reset();
    else
        store_.reset(new var_data<T,Str>(size()));
}

template <typename T, typename Str>
void var_acc<T,Str>::set_size(size_t size)
{
    current_ = bundle<T>(size, current_.target());
    if (valid())
        store_.reset(new var_data<T,Str>(size));
}

template <typename T, typename Str>
void var_acc<T,Str>::set_batch_size(uint64_t batch_size)
{
    if (!batch_size) throw std::invalid_argument("ALEA batch size must be positive");
    current_.target() = batch_size;
    current_.reset();
}

template <typename T, typename Str>
void var_acc<T,Str>::add(const computed<T> &source, uint64_t count,
                         var_acc<T,Str> *cascade)
{
    internal::check_valid(*this);
    if (source.size() != size()) throw size_mismatch();
    source.add_to(view<T>(current_.sum().data(), current_.size()));
    current_.count() += count;

    if (current_.is_full())
        add_bundle(cascade);
}

template <typename T, typename Str>
var_acc<T,Str> &var_acc<T,Str>::operator<<(const var_result<T,Str> &other)
{
    internal::check_valid(*this);
    if (size() != other.size())
        throw size_mismatch();

    // Leave partial bins in place and never mutate a const input result.
    auto incoming = other.store();
    incoming.unnormalize();
    internal::merge_moments(*store_, incoming, [](auto const& delta) {
        return delta.unaryExpr(typename bind<Str,T>::abs2_op()).eval();
    });
    return *this;
}

template <typename T, typename Str>
var_result<T,Str> var_acc<T,Str>::result() const
{
    internal::check_valid(*this);
    var_result<T,Str> result;
    var_acc<T,Str>(*this).finalize_to(result, nullptr);
    return result;
}

template <typename T, typename Str>
var_result<T,Str> var_acc<T,Str>::finalize()
{
    var_result<T,Str> result;
    finalize_to(result, nullptr);
    return result;
}

template <typename T, typename Str>
void var_acc<T,Str>::finalize_to(var_result<T,Str> &result, var_acc<T,Str> *cascade)
{
    internal::check_valid(*this);

    // add leftover data to the variance.  The upwards propagation must be
    // handled by going through a hierarchy in ascending order.
    if (current_.count() != 0)
        add_bundle(cascade);

    // data swap
    result.store_.reset();
    result.store_.swap(store_);

    // post-processing to result
    result.store_->normalize();
}

template <typename T, typename Str>
void var_acc<T,Str>::add_bundle(var_acc<T,Str> *cascade)
{
    internal::add_mean(*store_, (current_.sum() / current_.count()).eval(),
                       current_.count(), double(current_.count()) * current_.count(),
                       [](auto const& delta) { return delta.unaryExpr(typename bind<Str,T>::abs2_op()).eval(); });

    // add batch mean also to uplevel
    if (cascade != nullptr)
        cascade->add(make_adapter(current_.sum()), current_.count(), cascade+1);

    current_.reset();
}

template class var_acc<double>;
template class var_acc<std::complex<double>, circular_var>;
template class var_acc<std::complex<double>, elliptic_var>;

// We need an explicit copy constructor, as we need to copy the data
template <typename T, typename Str>
var_result<T,Str>::var_result(const var_result &other)
    : store_(other.store_ ? new var_data<T,Str>(*other.store_) : nullptr)
{ }

template <typename T, typename Str>
var_result<T,Str> &var_result<T,Str>::operator=(const var_result &other)
{
    store_.reset(other.store_ ? new var_data<T,Str>(*other.store_) : nullptr);
    return *this;
}

template <typename T, typename Strategy>
bool operator==(const var_result<T, Strategy> &r1, const var_result<T, Strategy> &r2)
{
    if (r1.count() == 0 && r2.count() == 0)
        return true;

    return r1.count() == r2.count()
        && r1.count2() == r2.count2()
        && r1.store().data() == r2.store().data()
        && r1.store().data2() == r2.store().data2();
}

template bool operator==(const var_result<double> &r1, const var_result<double> &r2);
template bool operator==(const var_result<std::complex<double>, circular_var> &r1,
                         const var_result<std::complex<double>, circular_var> &r2);
template bool operator==(const var_result<std::complex<double>, elliptic_var> &r1,
                         const var_result<std::complex<double>, elliptic_var> &r2);

template <typename T, typename Str>
column<typename var_result<T,Str>::var_type> var_result<T,Str>::stderror() const
{
    internal::check_valid(*this);
    return (store_->data2() / observations()).cwiseSqrt();
}

template <typename T, typename Str>
void var_result<T,Str>::reduce(const reducer &r)
{
    internal::check_reduction(r, valid(), valid() ? size() : 0);
    var_result staged(*this);
    staged.reduce_unchecked(r);
    store_.swap(staged.store_);
}

template <typename T, typename Str>
void var_result<T,Str>::reduce_unchecked(const reducer &r)
{
    internal::reduce_moments(*store_, r, [](auto const& delta) { return delta.unaryExpr(typename bind<Str,T>::abs2_op()).eval(); });
    if (!r.get_setup().have_result) store_.reset();
}

template class var_result<double>;
template class var_result<std::complex<double>, circular_var>;
template class var_result<std::complex<double>, elliptic_var>;


template <typename T, typename Str>
void serialize(serializer &s, const std::string &key, const var_result<T,Str> &self)
{
    internal::check_valid(self);
    internal::result_sentry group(s, key, internal::result_kind::variance);

    // serialize to uint64_t to make sure we are consistent across 32/64 bit
    serialize(s, "@size", static_cast<uint64_t>(self.store_->data_.size()));
    serialize(s, "count", self.store_->count_);
    serialize(s, "count2", self.store_->count2_);
    {
    internal::serializer_sentry subgroup(s, "mean");
    serialize(s, "value", self.store_->data_);
    serialize(s, "error", self.stderror());   // TODO temporary
    }
    serialize(s, "var", self.store_->data2_);
}

template <typename T, typename Str>
void deserialize(deserializer &s, const std::string &key, var_result<T,Str> &self)
{
    typedef typename var_result<T,Str>::var_type var_type;
    internal::result_reader_sentry group(s, key, internal::result_kind::variance);
    var_result<T,Str> staged;

    // deserialize from uint64_t
    uint64_t new_size_des;
    deserialize(s, "@size", new_size_des);
    if (!new_size_des || new_size_des > uint64_t(std::numeric_limits<Eigen::Index>::max())) throw size_mismatch();

    // first deserialize the fundamentals and make sure that the target fits
    {
        internal::deserializer_sentry mean(s, "mean");
        auto shape = s.get_shape("value");
        if (shape.size() != 1 || shape[0] != new_size_des) throw size_mismatch();
    }
    size_t new_size = new_size_des;
    if (!staged.valid() || staged.size() != new_size)
        staged.store_.reset(new var_data<T,Str>(new_size));

    // deserialize data
    deserialize(s, "count", staged.store_->count_);
    deserialize(s, "count2", staged.store_->count2_);
    if (!internal::valid_weight_count(staged.store_->count_, staged.store_->count2_))
        throw std::runtime_error("invalid ALEA squared-weight count");
    {
    internal::deserializer_sentry subgroup(s, "mean");
    deserialize(s, "value", staged.store_->data_);
    Eigen::Matrix<var_type, Eigen::Dynamic, 1> discard(staged.size());
    deserialize(s, "error", discard); // discard
    }
    deserialize(s, "var", staged.store_->data2_);

    self.store_.swap(staged.store_);
}

template void serialize(serializer &, const std::string &key, const var_result<double, circular_var> &);
template void serialize(serializer &, const std::string &key, const var_result<std::complex<double>, circular_var> &);
template void serialize(serializer &, const std::string &key, const var_result<std::complex<double>, elliptic_var> &);

template void deserialize(deserializer &, const std::string &key, var_result<double, circular_var> &);
template void deserialize(deserializer &, const std::string &key, var_result<std::complex<double>, circular_var> &);
template void deserialize(deserializer &, const std::string &key, var_result<std::complex<double>, elliptic_var> &);


template <typename T, typename Str>
std::ostream &operator<<(std::ostream &str, const var_result<T,Str> &self)
{
    internal::check_valid(self);
    internal::format_sentry sentry(str);
    verbosity verb = internal::get_format(str, PRINT_TERSE);

    if (verb == PRINT_VERBOSE)
        str << "<X> = ";
    str << self.mean() << " +- " << self.stderror();
    return str;
}

template std::ostream &operator<<(std::ostream &, const var_result<double, circular_var> &);
template std::ostream &operator<<(std::ostream &, const var_result<std::complex<double>, circular_var> &);
template std::ostream &operator<<(std::ostream &, const var_result<std::complex<double>, elliptic_var> &);

}}
