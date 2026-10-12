// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
// Independent cubic-polynomial interpolation and Fourier integral oracles for
// the retained Gull/Werner/Fuchs/Troyer spline. No external code adaptation.
#include "../fouriertransform.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <vector>

using dense_matrix=boost::numeric::ublas::matrix<double,boost::numeric::ublas::column_major>;
void generate_spline_matrix(dense_matrix&,double);
void evaluate_second_derivatives(double,dense_matrix&,std::vector<double>&,
                                std::vector<double>&,double,double,double);

namespace {
constexpr double beta=2., b=.3, c=.25, d=-.12;
constexpr double a=(-1.-b*beta-c*beta*beta-d*beta*beta*beta)/2.;
constexpr double c2=2*b+2*c*beta+3*d*beta*beta;
constexpr double c3=-(4*c+6*d*beta);
double cubic(double t) { return a+b*t+c*t*t+d*t*t*t; }
struct CubicTransformer: FourierTransformer {
  CubicTransformer():FourierTransformer(beta,1,1) {
    c2_[0][0][0]=c2; c3_[0][0][0]=c3;
  }
};
}

TEST(DmftSpline, TwoIntervalsRecoverExactCubicCurvatures) {
  const double dt=beta/2.;
  dense_matrix inverse(3,3);
  generate_spline_matrix(inverse,dt);
  std::vector<double> values{cubic(0),cubic(dt),cubic(beta)},curvature(3);
  evaluate_second_derivatives(dt,inverse,values,curvature,1.,c2,c3);
  for (int i=0;i<3;++i)
    EXPECT_NEAR(curvature[i],2*c+6*d*i*dt,1.e-13);
}

TEST(DmftSpline, CubicFourierIntegralMatchesAnalyticShortAndLargerMeshes) {
  const CubicTransformer transform;
  for (unsigned intervals : {2u,3u,8u}) {
    SCOPED_TRACE(intervals);
    itime_green_function_t time(intervals+1,1,1);
    for (unsigned i=0;i<=intervals;++i) time(i,0,0,0)=cubic(beta*i/intervals);
    matsubara_green_function_t frequency(7,1,1);
    transform.forward_ft(time,frequency);
    for (unsigned n=0;n<frequency.nfreq();++n) {
      const std::complex<double> iw(0,std::acos(-1.)*(2*n+1)/beta);
      // Integrating the cubic against exp(iw*t) by parts terminates after
      // its third derivative: endpoint sums fix all four exact coefficients.
      const auto exact=1./iw+c2/(iw*iw)+c3/(iw*iw*iw)+12*d/(iw*iw*iw*iw);
      EXPECT_NEAR(std::abs(frequency(n,0,0,0)-exact),0.,1.e-12);
    }
  }
}
