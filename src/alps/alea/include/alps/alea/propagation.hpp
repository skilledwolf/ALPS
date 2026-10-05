/*
 * Copyright (C) 1998-2018 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 * For use in publications, see ACKNOWLEDGE.TXT
 */
#pragma once

#include <Eigen/Core>
#include <cmath>
#include <stdexcept>

#include <alps/alea/complex_op.hpp>

// TODO maybe a better way?
#include <alps/alea/mean.hpp>
#include <alps/alea/variance.hpp>
#include <alps/alea/covariance.hpp>
#include <alps/alea/autocorr.hpp>
#include <alps/alea/batch.hpp>

namespace alps { namespace alea {

/**
 * Do not perform error propagation.
 *
 * Given a transformation `f` and random sample `X`, just trasform the sample
 * mean `f(Mean[X])` and discard any error information.
 */
struct no_prop { };

/**
 * Perform linearized error propagation by estimating the Jacobian.
 *
 * Given a transformation `f` and a random sample `X`, estimate the propagated
 * uncertainties by performing a Taylor series and keeping the linear term:
 *
 *     Cov[f(X)] = df/dX Cov[X] (df/dX)^T + O(d^2f/dx^2)
 *
 * where `df/dX` is the Jacobian of `f` at `X`, as estimated by finite
 * central differences of `dx`. A zero step chooses a relative step for each
 * component from machine precision, independently of the sample's error.
 * This procedure is exact for declared linear transformations;
 * for non-linear transformation, it will introduce bias.
 *
 * @see alps::alea::jacobian
 */
struct linear_prop
{
    linear_prop() : dx_(0) { }

    linear_prop(double dx) : dx_(dx) {
        if (!std::isfinite(dx) || dx < 0)
            throw std::invalid_argument("Jacobian step must be finite and nonnegative");
    }

    double dx() const { return dx_; }

private:
    double dx_;
};

/**
 * Perform Jackknife rebatching.
 *
 * Jackknife is a rebatching method, which can operate on any distribution and
 * exactly removes the bias in the transformed uncertainties up to order `1/N`,
 * where `N` is the sample size.
 *
 * @see alps::alea::jackknife
 */
struct jackknife_prop { };

/**
 * Given a function `f`, estimate its Jacobian `J[i,j] = df[i]/dx[j]`.
 *
 * Estimate the Jacobian of a transformation `f` at the point `x` by central
 * differences:
 *
 *       J[i,j] ~= (f(x + dx e[j]) - f(x - dx e[j]))[i] / (2 dx);
 *
 * where `e[j]` denotes the `j`-th unit vector. A zero `dx` chooses a relative
 * step per component; declared linear transformers are evaluated on the
 * basis vectors instead. Transformer exceptions propagate unchanged.
 */
template <typename T>
typename eigen<T>::matrix jacobian(const transformer<T> &f, column<T> x, double dx);


/**
 * Perform Jackknife transformation to weighted pseudovalues. Empty bins are
 * skipped, and at least two occupied bins are required. Transformer exceptions
 * propagate unchanged, including singular leave-one-out estimates.
 */
template <typename T>
batch_data<T> jackknife(const batch_data<T> &in, const transformer<T> &tf);

}}
