/*
 * Copyright (C) 1998-2018 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 * For use in publications, see ACKNOWLEDGE.TXT
 */

#pragma once

#include <alps/alea/core.hpp>
#include <alps/alea/util.hpp>
#include <alps/alea/variance.hpp>
#include <alps/alea/covariance.hpp>
#include <alps/alea/fisher_f.hpp>
#include <alps/alea/internal/util.hpp>

#include <algorithm>
#include <cmath>

namespace alps { namespace alea {
    class t2_result;
}}

namespace alps { namespace alea {

/**
 * Result of Hotelling's T^2 test
 *
 * @see test_mean(), t2_test()
 */
class t2_result
{
public:
    typedef fisher_f_distribution dist_type;

public:
    /** Initializes t2 test */
    t2_result(double score, double f_size, double f_dof)
        : score_(score), dist_(f_size, f_dof)
    { }

    /** Returns whether the lower alternate hypothesis can be tested.
     *
     * For low values of size, the F distribution does not have a mode,
     * because there is no way to distinguish the lower alternate hypothesis
     * (too large error bars) from an accidental "hit".
     */
    bool has_plower() const { return dist_.degrees_of_freedom1() > 2; }

    /** p-value in favour of the lower alternate Hypothesis */
    double pvalue_lower() const { return has_plower() ? dist_.cdf(score_) : 1; }

    /** p-value in favour of the upper alternate Hypothesis */
    double pvalue_upper() const { return dist_.ccdf(score_); }

    /** Standard upper-tail p-value for rejecting equality of means. */
    double pvalue() const { return pvalue_upper(); }

    /** t2 statistic and argument to `dist()` */
    double score() const { return score_; }

    /** distribution */
    const dist_type &dist() const { return dist_; }

private:
    double score_;
    dist_type dist_;
};

/**
 * Perform Hotelling's T^2 test given a set of uncorrelated differences.
 *
 * `var` is the variance of each mean difference, after diagonalizing its
 * covariance. `covariance_dof` is the degrees of freedom of the covariance
 * estimate: n-1 for one sample, n1+n2-2 for two equal-covariance samples.
 * `atol` is the absolute tolerance for deterministic directions.
 *
 * @see test_mean()
 */
t2_result t2_test(const column<double> &diff, const column<double> &var,
                  double covariance_dof, double atol=1e-14);

namespace internal {
t2_result test_covariance(const column<double>& diff, const Eigen::MatrixXd& cov,
                          double dof, double atol);

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

template<class Result> double test_observations(const Result& result)
{
    if (!result.valid() || result.count() == 0)
        throw std::invalid_argument("Mean test requires a nonempty result");
    if (!valid_weight_count(result.count(), result.count2()))
        throw std::invalid_argument("Mean test requires consistent observation weights");
    double n = result.observations();
    if (!std::isfinite(n) || n <= 1 || n > double(result.count()))
        throw std::invalid_argument("Mean test requires more than one independent observation");
    return n;
}
}

/**
 * Test a Gaussian sample mean against a known vector. Effective observation
 * counts give an approximation for weighted bins or correlated samples.
 * Complex results are tested as joint real/imaginary vectors; raw batches or
 * elliptic covariance are required, without a circularity assumption.
 */
template <typename Result, typename Derived>
t2_result test_mean(const Result &result,
                    const Eigen::MatrixBase<Derived> &expected,
                    double atol=1e-14)
{
    static_assert(is_alea_result<Result>::value, "Result is not alea result");
    static_assert(traits<Result>::HAVE_VAR, "Result1 must have variance");

    double n = internal::test_observations(result);
    auto diff = internal::real_vector(result.mean());
    auto target = internal::real_vector(expected);
    if (diff.size() != target.size()) throw size_mismatch();
    diff -= target;
    return internal::test_covariance(diff, internal::real_covariance(result)/n, n-1, atol);
}

/**
 * Test mean of stochastic result `result` against known result `expected`.
 */
template <typename Result, typename Derived>
t2_result test_mean(const Eigen::MatrixBase<Derived> &expected,
                    const Result &result,
                    double atol=1e-14)
{
    return test_mean(result, expected, atol);
}

/**
 * Test two independent Gaussian sample means with equal population covariance.
 * Their sample covariances are pooled before evaluating Hotelling's statistic.
 * Effective observation counts give an approximation for weighted bins.
 */
template <typename Result1, typename Result2>
typename std::enable_if<
    is_alea_result<Result1>::value && is_alea_result<Result2>::value,
    t2_result>::type
test_mean(const Result1 &result1, const Result2 &result2, double atol=1e-14)
{
    static_assert(traits<Result1>::HAVE_VAR, "Result1 must have variance");
    static_assert(traits<Result2>::HAVE_VAR, "Result2 must have variance");

    double n1 = internal::test_observations(result1), n2 = internal::test_observations(result2);
    auto diff = internal::real_vector(result1.mean());
    auto other = internal::real_vector(result2.mean());
    if (diff.size() != other.size()) throw size_mismatch();
    diff -= other;
    double dof = n1+n2-2;
    Eigen::MatrixXd covariance = ((n1-1)*internal::real_covariance(result1)
                               + (n2-1)*internal::real_covariance(result2))/dof;
    covariance *= 1/n1 + 1/n2;
    return internal::test_covariance(diff, covariance, dof, atol);
}


}} /* namespace alps::alea */
