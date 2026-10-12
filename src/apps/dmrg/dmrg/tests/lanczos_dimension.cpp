// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
// Independent diagonal operators expose an N-1 Krylov truncation: the seed has
// overlap with every distinct eigenvalue, so exact recovery can require N steps.
#define WITH_LAPACK
#include <valarray>
#include "dmtk/lanczos.cc"
#include <alps/testing/temporary_directory.hpp>
#include <gtest/gtest.h>
#include <complex>

namespace reference {
struct diagonal { std::vector<double> values; };
template<class T> dmtk::Vector<T> product(const diagonal& matrix,const dmtk::Vector<T>& vector) {
    dmtk::Vector<T> output(vector.size());
    for(std::size_t i=0;i<vector.size();++i) output[i]=matrix.values[i]*vector[i];
    return output;
}
template<class T> void complete_space() {
    alps::testing::TemporaryDirectory tmp;
    for(const int dimension:{2,3,6}) {
        SCOPED_TRACE(dimension);
        diagonal matrix; for(int i=0;i<dimension;++i) matrix.values.push_back(double(i));
        dmtk::Vector<T> seed(dimension); seed=T(1);
        dmtk::Vector<dmtk::Vector<T>> states(1); states[0].resize(dimension);
        dmtk::Vector<double> energies(1),a(100),b(100); int iterations=-1;
        const auto filename=(tmp.path()/"vectors.dat").string();
        dmtk::lanczos<T,diagonal,dmtk::Vector<T>>(matrix,states,seed,energies,1,a,b,
            iterations,1.e-14,true,true,filename.c_str());
        EXPECT_EQ(iterations,dimension); EXPECT_NEAR(energies[0],0.,1.e-11);
        double norm=0,residual=0;
        for(int i=0;i<dimension;++i) {
            norm+=std::norm(states[0][i]);
            residual+=std::norm((matrix.values[i]-energies[0])*states[0][i]);
        }
        EXPECT_NEAR(norm,1.,1.e-11); EXPECT_LT(std::sqrt(residual),1.e-10);
    }
}
}
TEST(DMTKLanczos, RealTinySectorsUseCompleteKrylovSpace) { reference::complete_space<double>(); }
TEST(DMTKLanczos, ComplexTinySectorsUseCompleteKrylovSpace) { reference::complete_space<std::complex<double>>(); }
