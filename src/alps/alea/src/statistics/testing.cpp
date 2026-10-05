/*
 * Copyright (C) 1998-2018 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 * For use in publications, see ACKNOWLEDGE.TXT
 */
#include <alps/alea/testing.hpp>
#include <Eigen/Eigenvalues>
#include <limits>

namespace alps { namespace alea {
t2_result t2_test(const column<double>& diff, const column<double>& var,
                  double dof, double atol)
{
    if (diff.size() != var.size()) throw size_mismatch();
    if (diff.size() == 0) throw std::invalid_argument("Mean test requires a nonempty vector");
    if (!std::isfinite(dof) || dof <= 0 || !std::isfinite(atol) || atol < 0 ||
        !diff.allFinite() || !var.allFinite())
        throw std::invalid_argument("Invalid mean-test data, covariance degrees of freedom or tolerance");
    // Sample-induced singularity cannot be treated as a deterministic
    // population constraint. Check the full dimension before reducing rank.
    if (dof <= diff.size()-1)
        throw std::invalid_argument("Too few independent observations for the tested vector dimension");
    if ((var.array() < 0).any())
        throw std::invalid_argument("Mean-test covariance must be positive semidefinite");
    double t2 = 0;
    size_t rank = 0;
    bool deterministic_mismatch = false;
    for (Eigen::Index i=0; i<diff.size(); ++i) {
        if (var(i) == 0) {
            deterministic_mismatch |= std::abs(diff(i)) > atol;
        } else {
            double standardized = diff(i)/std::sqrt(var(i));
            t2 += standardized*standardized;
            ++rank;
        }
    }
    if (!rank) return t2_result(deterministic_mismatch ? INFINITY : 0., 1, dof);
    double denominator_dof = dof-rank+1;
    double score = deterministic_mismatch ? INFINITY : denominator_dof/(rank*dof)*t2;
    return t2_result(score, rank, denominator_dof);
}

namespace internal {
t2_result test_covariance(const column<double>& diff, const Eigen::MatrixXd& cov,
                          double dof, double atol)
{
    if (cov.rows() != diff.size() || cov.cols() != diff.size()) throw size_mismatch();
    if (!cov.allFinite() || !cov.isApprox(cov.transpose(), 1e-12))
        throw std::invalid_argument("Mean-test covariance must be finite and symmetric");
    if (diff.size() == 0) throw std::invalid_argument("Mean test requires a nonempty vector");
    // Diagonalize a correlation matrix so component units cannot change the
    // numerical rank. Divide rows/columns separately to avoid scale overflow.
    column<double> scale = cov.diagonal().cwiseSqrt();
    Eigen::MatrixXd correlation = cov;
    for (Eigen::Index i=0; i<scale.size(); ++i) {
        if (cov(i,i) < 0) throw std::invalid_argument("Mean-test covariance has negative variance");
        if (scale(i) == 0) {
            if (!cov.row(i).isZero(0)) throw std::invalid_argument("Zero variance has nonzero covariance");
            scale(i) = 1;
        }
        correlation.row(i) /= scale(i);
        correlation.col(i) /= scale(i);
    }
    if (!correlation.allFinite() || !correlation.isApprox(correlation.transpose(), 1e-12))
        throw std::invalid_argument("Mean-test correlation must be finite and symmetric");
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen(correlation);
    if (eigen.info() != Eigen::Success)
        throw std::invalid_argument("Mean-test covariance diagonalization failed");
    column<double> variance = eigen.eigenvalues();
    double tolerance = variance.cwiseAbs().maxCoeff()*std::numeric_limits<double>::epsilon()*variance.size();
    for (Eigen::Index i=0; i<variance.size(); ++i)
        if (std::abs(variance(i)) <= tolerance) variance(i) = 0;
    return t2_test(eigen.eigenvectors().transpose()*diff.cwiseQuotient(scale), variance, dof, atol);
}
}
}} /* namespace alps::alea */
