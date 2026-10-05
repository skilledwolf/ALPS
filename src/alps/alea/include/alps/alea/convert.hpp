/*
 * Copyright (C) 1998-2018 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 * For use in publications, see ACKNOWLEDGE.TXT
 */
#pragma once

#include <alps/alea/core.hpp>
#include <alps/alea/mean.hpp>
#include <alps/alea/covariance.hpp>
#include <alps/alea/autocorr.hpp>
#include <alps/alea/batch.hpp>

#include <alps/alea/internal/util.hpp>
#include <alps/alea/internal/joined.hpp>

// Forward

namespace alps { namespace alea {
    template <typename Result> struct joiner;

    template <typename T> class mean_dat;
    template <typename T, typename Str> class var_data;
    template <typename T, typename Str> class cov_data;
    template <typename T> class batch_data;
}}

// Actual

namespace alps { namespace alea {

/**
 * Joins two statistical results together, assuming their mutual independence.
 *
 * Returns a result of the random vector which is the concatenation of the
 * random vectors corresponding to the arguments `first`, `second`, i.e., of
 * size ` first.size() + second.size()`.   Assumes that `first` and `second`
 * are uncorrelated (their covariance is zero).
 *
 * One can combine results of different type; in this case, `Result` is
 * inferred in such a way that as much information as possible is preserved
 * from the constituent accumulators.
 *
 * @see alps::alea::internal::joined
 */
template <typename R1, typename R2,
          typename Result=typename internal::joined<R1, R2>::result_type>
Result join(const R1 &first, const R2 &second)
{
    return joiner<Result>()(first, second);
}

/** Helper class for joining results together */
template <typename Result>
struct joiner;

template <typename T>
struct joiner<mean_result<T> >
{
    template <typename R1, typename R2>
    mean_result<T> operator()(const R1 &first, const R2 &second)
    {
        if (first.store().count() != second.store().count())
            throw weight_mismatch();  // TODO

        mean_result<T> res(mean_data<T>(first.size() + second.size()));
        res.store().data().topRows(first.size()) = first.store().data();
        res.store().data().bottomRows(second.size()) = second.store().data();
        res.store().count() = first.store().count();
        return res;
    }
};

template <typename T, typename Str>
struct joiner<var_result<T,Str> >
{
    template <typename R1, typename R2>
    var_result<T,Str> operator()(const R1 &first, const R2 &second)
    {
        if (first.store().count() != second.store().count())
            throw weight_mismatch();
        if (first.store().count2() != second.store().count2())
            throw weight_mismatch();

        var_result<T,Str> res(var_data<T,Str>(first.size() + second.size()));
        res.store().data().topRows(first.size()) = first.store().data();
        res.store().data().bottomRows(second.size()) = second.store().data();
        res.store().data2().topRows(first.size()) = first.store().data2();
        res.store().data2().bottomRows(second.size()) = second.store().data2();
        res.store().count() = first.store().count();
        res.store().count2() = first.store().count2();
        return res;
    }
};

template <typename T, typename Str>
struct joiner<cov_result<T,Str> >
{
    template <typename R1, typename R2>
    cov_result<T,Str> operator()(const R1 &first, const R2 &second)
    {
        if (first.store().count() != second.store().count())
            throw weight_mismatch();
        if (first.store().count2() != second.store().count2())
            throw weight_mismatch();

        cov_result<T,Str> res(cov_data<T,Str>(first.size() + second.size()));
        res.store().data().topRows(first.size()) = first.store().data();
        res.store().data().bottomRows(second.size()) = second.store().data();

        // ignore cross correlation
        res.store().data2().topLeftCorner(first.size(), first.size())
                                                = first.store().data2();
        res.store().data2().bottomRightCorner(second.size(), second.size())
                                                = second.store().data2();
        res.store().count() = first.store().count();
        res.store().count2() = first.store().count2();
        return res;
    }
};

template <typename T>
struct joiner<autocorr_result<T> >
{
    template <typename R1, typename R2>
    autocorr_result<T> operator()(const R1 &first, const R2 &second)
    {
        if (first.count() != second.count())
            throw weight_mismatch();
        if (first.nlevel() != second.nlevel())
            throw size_mismatch();

        // granularities are checked on the individual levels
        autocorr_result<T> res(first.nlevel());
        for (size_t l = 0; l != first.nlevel(); ++l)
            res.level(l) = join(first.level(l), second.level(l));
        return res;
    }
};

template <typename T>
struct joiner<batch_result<T> >
{
    template <typename R1, typename R2>
    batch_result<T> operator()(const R1 &first, const R2 &second)
    {
        if (first.store().count() != second.store().count())
            throw weight_mismatch();
        if (first.store().num_batches() != second.store().num_batches())
            throw size_mismatch();

        batch_result<T> res(batch_data<T>(first.size() + second.size(),
                                          first.store().num_batches()));

        res.store().batch().topRows(first.size()) = first.store().batch();
        res.store().batch().bottomRows(second.size()) = second.store().batch();
        res.store().count() = first.store().count();
        return res;
    }
};

namespace internal {
template<class Derived>
column<double> real_vector(const Eigen::MatrixBase<Derived>& value)
{
    if (value.rows() != 1 && value.cols() != 1) throw size_mismatch();
    constexpr bool complex = Eigen::NumTraits<typename Derived::Scalar>::IsComplex;
    column<double> out(value.size() * (complex ? 2 : 1));
    for (Eigen::Index i=0; i<value.size(); ++i) {
        if constexpr (complex) {
            out(2*i) = value.derived().coeff(i).real();
            out(2*i+1) = value.derived().coeff(i).imag();
        } else out(i) = value.derived().coeff(i);
    }
    return out;
}

template<class Result>
Eigen::MatrixXd real_covariance(const Result& result)
{
    if constexpr (!Eigen::NumTraits<typename traits<Result>::value_type>::IsComplex) {
        if constexpr (traits<Result>::HAVE_COV) return result.cov();
        else {
            if (result.size() != 1)
                throw std::invalid_argument("Multicomponent mean tests require full covariance");
            return result.var().asDiagonal();
        }
    } else {
        // Circular covariance does not determine real/imaginary covariance.
        // Raw batches and elliptic results do retain the necessary information.
        if constexpr (traits<Result>::HAVE_BATCH ||
                      std::is_same_v<typename traits<Result>::var_type, complex_op<double>>) {
            auto blocks = [&] {
                if constexpr (traits<Result>::HAVE_BATCH)
                    return result.template cov<elliptic_var>();
                else if constexpr (traits<Result>::HAVE_COV) return result.cov();
                else {
                    if (result.size() != 1)
                        throw std::invalid_argument("Multicomponent complex mean tests require full covariance");
                    return typename eigen<complex_op<double>>::matrix(result.var().asDiagonal());
                }
            }();
            Eigen::MatrixXd out(2*result.size(), 2*result.size());
            for (size_t i=0; i<result.size(); ++i) for (size_t j=0; j<result.size(); ++j) {
                auto b = blocks(i,j);
                out(2*i,2*j)=b.rere(); out(2*i,2*j+1)=b.reim();
                out(2*i+1,2*j)=b.imre(); out(2*i+1,2*j+1)=b.imim();
            }
            return out;
        } else throw std::invalid_argument("Complex mean tests require raw batches or elliptic covariance");
    }
}

}

/** Expose retained real/imaginary evidence as interleaved real components.
 * Circular complex summary moments cannot reconstruct elliptic covariance.
 */
template<class R> auto real_components(R const& in) {
    static_assert(Eigen::NumTraits<typename traits<R>::value_type>::IsComplex,
                  "real_components expects a complex result");
    internal::check_valid(in);
    if constexpr (traits<R>::HAVE_BATCH) {
        batch_data<double> data(2*in.size(),in.num_batches());
        for (size_t i=0;i<in.num_batches();++i)
            data.batch().col(i)=internal::real_vector(in.store().batch().col(i));
        data.count()=in.store().count();
        return batch_result<double>(data);
    } else if constexpr (traits<R>::HAVE_VAR) {
        cov_data<double> data(2*in.size());
        data.data()=internal::real_vector(in.mean());
        data.data2()=internal::real_covariance(in);
        data.count()=in.count(); data.count2()=in.count2();
        return cov_result<double>(data);
    } else {
        mean_data<double> data(2*in.size());
        data.data()=internal::real_vector(in.mean()); data.count()=in.count();
        return mean_result<double>(data);
    }
}

/** Pool independent runs without inventing a resumable combined time series. */
template<class R> R merge(std::vector<R> const& runs) {
    if (runs.empty()) throw std::invalid_argument("ALEA merge requires at least one result");
    using T=typename traits<R>::value_type;
    auto size=runs.front().valid() ? runs.front().size() : 0;
    uint64_t count=0;
    for (auto const& run:runs) {
        internal::check_valid(run);
        if (!size || run.size()!=size) throw size_mismatch();
        if (run.count()>std::numeric_limits<uint64_t>::max()-count)
            throw std::overflow_error("ALEA merged sample count overflows");
        count+=run.count();
    }
    if constexpr (traits<R>::HAVE_BATCH) {
        size_t slots=0;
        auto limit=size_t(std::numeric_limits<Eigen::Index>::max())/size;
        for (auto const& run:runs) {
            if (run.num_batches()>limit-slots) throw size_mismatch();
            slots+=run.num_batches();
        }
        batch_data<T> data(size,slots);
        size_t offset=0;
        for (auto const& run:runs) {
            data.batch().middleCols(offset,run.num_batches())=run.store().batch();
            data.count().segment(offset,run.num_batches())=run.store().count();
            offset+=run.num_batches();
        }
        return R(data);
    } else if constexpr (traits<R>::HAVE_TAU) {
        size_t depth=std::numeric_limits<size_t>::max();
        for (auto const& run:runs) if (run.count()) depth=std::min(depth,run.nlevel());
        if (!count) depth=1;
        R result(depth);
        for (size_t i=0;i<depth;++i) {
            std::vector<var_result<T>> levels;
            for (auto const& run:runs) if (run.count()) levels.push_back(run.level(i));
            result.level(i)=count ? merge(levels) : var_acc<T>(size).result();
        }
        return result;
    } else if constexpr (traits<R>::HAVE_COV) {
        cov_acc<T,typename traits<R>::strategy_type> acc(size);
        for (auto const& run:runs) acc << run;
        return acc.result();
    } else if constexpr (traits<R>::HAVE_VAR) {
        var_acc<T,typename traits<R>::strategy_type> acc(size);
        for (auto const& run:runs) acc << run;
        return acc.result();
    } else {
        mean_acc<T> acc(size);
        for (auto const& run:runs) acc << run;
        return acc.result();
    }
}

}}
