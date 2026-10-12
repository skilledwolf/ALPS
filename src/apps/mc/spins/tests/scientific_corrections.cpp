// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
// Independent Hamiltonian evidence for corrections to original ALPS
// Troyer/Stoeckli spinmc matrix parsing and onsite-energy measurements.
#include "factory.h"
#include "spinsim.h"
#include <memory>
#include "ising.h"
#include "on.h"
#include <boost/math/special_functions/bessel.hpp>
#include <gtest/gtest.h>
#include <iomanip>
#include <sstream>

namespace {
alps::Parameters parameters(const std::string &model) {
    alps::Parameters result;
    result["MODEL"] = model; result["LATTICE"] = "dimer";
    result["LATTICE_LIBRARY"] = std::string(ALPS_SPINMC_LATTICES);
    result["UPDATE"] = "local"; result["T"] = 1.2;
    result["SWEEPS"] = 100000; result["THERMALIZATION"] = 0;
    result["SEED"] = 173; result["RNG"] = "mt19937";
    result["J"] = "0"; result["h"] = "0"; result["S"] = 1.5;
    return result;
}
struct original_run {
    explicit original_run(const alps::Parameters &input)
        : worker(SpinFactory().make_worker(alps::ProcessList(1), input, 0)),
          measurements(dynamic_cast<alps::scheduler::MCRun *>(worker.get())) {
        if (!measurements) throw std::logic_error("missing original spinmc observation boundary");
        worker->start_worker();
        // The scheduler clears measurements on its first thermalized boundary.
        if (int(input["THERMALIZATION"]) == 0) worker->run();
    }
    std::unique_ptr<alps::scheduler::Worker> worker;
    alps::scheduler::MCRun *measurements;
};
}

TEST(SpinMCScientificCorrections, ScalarAndDiagonalSquareInputsRepresentDeclaredHamiltonian) {
    alps::Parameters symbols; symbols["d"] = .7;
    std::string scalar = "d";
    SquareMatrix<double, 2> identity(scalar, symbols);
    EXPECT_DOUBLE_EQ(identity.getElement(0, 0), .7); EXPECT_DOUBLE_EQ(identity.getElement(1, 1), .7);
    EXPECT_DOUBLE_EQ(identity.getElement(0, 1), 0); EXPECT_DOUBLE_EQ(identity.getElement(1, 0), 0);
    EXPECT_NEAR(identity.vec_mat_vec(TinyVector<double, 2>(.6, .8)), .7, 1.e-15);
    std::string diagonal = "2 3 5";
    SquareMatrix<double, 3> anisotropy(diagonal, symbols);
    EXPECT_DOUBLE_EQ(anisotropy.det(), 30);
    // v=(1,2,3) gives vᵀDv=2+12+45, independent of parser indexing.
    EXPECT_DOUBLE_EQ(anisotropy.vec_mat_vec(TinyVector<double, 3>(1, 2, 3)), 59);
    for (unsigned i = 0; i < 3; ++i) for (unsigned j = 0; j < 3; ++j)
        EXPECT_DOUBLE_EQ(anisotropy.getElement_nc(i, j), i != j ? 0 : i == 0 ? 2 : i == 1 ? 3 : 5);
    std::string full4 = "1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16";
    SquareMatrix<double, 4> four(full4, symbols);
    TinyVector<double, 4> x, y;
    for (unsigned i = 0; i < 4; ++i) x[i] = i + 1;
    y[0] = -1; y[1] = .5; y[2] = 2; y[3] = -2;
    // Literal rows give m*y=(-2,-4,-6,-8), m*x=(30,70,110,150).
    EXPECT_DOUBLE_EQ(four.vec_mat_vec(x, y), -60);
    EXPECT_DOUBLE_EQ(four.vec_mat_vec(x), 1100);
}

TEST(SpinMCScientificCorrections, FullAndTriangularInputsRetainOffDiagonalHamiltonian) {
    alps::Parameters symbols;
    std::string triangular = "2 .3 -1", complete = "2 .3 .3 -1";
    SquareMatrix<double, 2> first(triangular, symbols), second(complete, symbols);
    ONMoment<2> spin; spin.state_ = TinyVector<double, 2>(.6, .8);
    const TinyVector<double, 2> reflection(1, 0);
    // H=-[2*x²+0.6*x*y-y²], at (+.6,.8) and (-.6,.8).
    EXPECT_NEAR(onsite_energy(spin, first), -.368, 1.e-15);
    EXPECT_NEAR(onsite_energy_change(spin, first, reflection), .576, 1.e-15);
    EXPECT_NEAR(onsite_energy(spin, second), -.368, 1.e-15);
    auto reflected = spin; reflected.update(reflection);
    EXPECT_NEAR(onsite_energy(reflected, second), .208, 1.e-15);
    EXPECT_NEAR(onsite_energy(reflected, second) - onsite_energy(spin, second), .576, 1.e-15);
    for (unsigned i = 0; i < 2; ++i) for (unsigned j = 0; j < 2; ++j)
        EXPECT_DOUBLE_EQ(first.getElement_nc(i, j), second.getElement(i, j));
}

TEST(SpinMCScientificCorrections, HugeFixedNormConstantsDoNotEnterVariableEnergyOrProposalActions) {
    alps::Parameters symbols;
    ONMoment<2> spin;spin.state_=TinyVector<double,2>(.6,.8+2.e-16);
    const TinyVector<double,2> direction(1.+2.e-16,0.);
    for(const char* offset:{"1e10","1e100"}) {
        std::string isotropic=std::string(offset)+" 0 0 "+offset;
        SquareMatrix<double,2> full(isotropic,symbols);
        const auto parts=onsite_energy_parts(spin,full);
        EXPECT_DOUBLE_EQ(parts.first,0.);
        EXPECT_DOUBLE_EQ(parts.second,-std::stod(offset));
        EXPECT_DOUBLE_EQ(onsite_energy(spin,full),parts.second);
        EXPECT_DOUBLE_EQ(onsite_energy_change(spin,full,direction),0.);
        EXPECT_DOUBLE_EQ(self_bond_energy_parts(spin,full).first,0.);
        MIdMatrix<double,2> scalar(std::stod(offset));
        const auto scalar_site=onsite_energy_parts(spin,scalar);
        const auto scalar_self=self_bond_energy_parts(spin,scalar);
        EXPECT_DOUBLE_EQ(scalar_site.first,0.);
        EXPECT_DOUBLE_EQ(scalar_site.second,-std::stod(offset));
        EXPECT_DOUBLE_EQ(scalar_self.first,0.);
        EXPECT_DOUBLE_EQ(scalar_self.second,-std::stod(offset));
    }
    std::string anisotropic="10000000000 .3 .3 10000000001";
    SquareMatrix<double,2> matrix(anisotropic,symbols);
    spin.state_=TinyVector<double,2>(.6,.8);
    // H=-1e10 - y² - .6*x*y: variable energy -.64-.288.
    EXPECT_NEAR(onsite_energy_parts(spin,matrix).first,-.928,2.e-15);
    EXPECT_NEAR(onsite_energy_change(spin,matrix,TinyVector<double,2>(1,0)),.576,2.e-15);
    std::string antisymmetric="0 1e100 -1e100 0";
    SquareMatrix<double,2> self(antisymmetric,symbols);
    EXPECT_DOUBLE_EQ(self_bond_energy_parts(spin,self).first,0.);
    EXPECT_DOUBLE_EQ(self_bond_energy_parts(spin,self).second,0.);
}

TEST(SpinMCScientificCorrections, IsotropicONOnsiteEnergyDependsOnDeclaredScalar) {
    for (const auto coefficient : {0., .7, -2.}) {
        MIdMatrix<double, 2> xy(coefficient); MIdMatrix<double, 3> heisenberg(coefficient);
        MIdMatrix<double, 4> o4(coefficient);
        EXPECT_DOUBLE_EQ(onsite_energy(ONMoment<2>{}, xy), -coefficient);
        EXPECT_DOUBLE_EQ(onsite_energy(ONMoment<3>{}, heisenberg), -coefficient);
        EXPECT_DOUBLE_EQ(onsite_energy(ONMoment<4>{}, o4), -coefficient);
    }
}

TEST(SpinMCScientificCorrections, EverySiteContributesItsConstantOnsiteEnergyInOriginalWorker) {
    for (const auto *model : {"Ising", "XY", "Heisenberg", "O(4)"}) {
        SCOPED_TRACE(model);
        auto input = parameters(model); input["D"] = ".7";
        original_run run(input);
        // Original Ising convention has +D; original O(N) convention has
        // -sᵀDs. Both include all 2 sites and the declared S²=2.25.
        const double exact = std::string(model) == "Ising" ? 3.15 : -3.15;
        for (unsigned update = 0; update < 128; ++update) {
            run.worker->run();
            EXPECT_NEAR(alps::RealObsevaluator(run.measurements->get_measurements()["Energy"]).mean(), exact, 1.e-11);
            EXPECT_NEAR(alps::RealObsevaluator(run.measurements->get_measurements()["Energy^2"]).mean(), exact * exact, 1.e-11);
        }
        EXPECT_NEAR(alps::RealObsevaluator(run.measurements->get_measurements()["Energy"]).mean(), exact, 1.e-11);
    }
}

TEST(SpinMCScientificCorrections, ForcedSquareStoragePreservesScalarAndDiagonalOnsiteHamiltonians) {
    for (const auto *model : {"XY", "Heisenberg"}) for (bool diagonal_input : {false, true}) {
        SCOPED_TRACE(std::string(model) + (diagonal_input ? " diagonal" : " scalar"));
        auto input = parameters(model);
        // Full zero J forces general matrix storage, independently of the
        // scalar/diagonal representation used for the same D=.7I.
        input["J"] = std::string(model) == "XY" ? "0 0 0 0" : "0 0 0 0 0 0 0 0 0";
        input["D"] = diagonal_input ? (std::string(model) == "XY" ? ".7 .7" : ".7 .7 .7") : ".7";
        original_run run(input);
        for (unsigned update = 0; update < 128; ++update) {
            run.worker->run();
            EXPECT_NEAR(alps::RealObsevaluator(run.measurements->get_measurements()["Energy"]).mean(), -3.15, 1.e-11);
        }
    }
}

TEST(SpinMCScientificCorrections, ScalarCouplingInForcedSquareStorageMatchesIndependentXYDimerPartition) {
    auto input = parameters("XY"); input["S"] = 1.; input["J"] = ".6";
    input["D"] = "0 0 0"; input["THERMALIZATION"] = 1000;
    original_run run(input);
    for (unsigned warmup=0; warmup<1000; ++warmup) run.worker->run();
    for (unsigned update = 0; update < 100000; ++update) {
        run.worker->run();
    }
    const double correlation = boost::math::cyl_bessel_i(1, .5) / boost::math::cyl_bessel_i(0, .5);
    EXPECT_NEAR(alps::RealObsevaluator(run.measurements->get_measurements()["Energy"]).mean(), -.6 * correlation, .025);
    EXPECT_NEAR(alps::RealObsevaluator(run.measurements->get_measurements()["Energy^2"]).mean(), .36 * (1 - 2 * correlation), .025);
}

TEST(SpinMCScientificCorrections, NonsymmetricExchangeProposalUsesOriginalBondOrientationAtBothEndpoints) {
    const auto input = parameters("XY");
    const alps::graph_helper<> graph(input);
    std::string text = "1 .7 -.4 -.2";
    SquareMatrix<double,2> matrix(text,input);
    const auto sources = spinmc_local_bond_sources<SquareMatrix<double,2>>(graph.graph());
    const auto bond = *graph.bonds().first;
    const auto left = graph.source(bond), right = graph.target(bond);
    ONMoment<2> x,y;
    x.state_ = TinyVector<double,2>(.6,.8);
    y.state_ = TinyVector<double,2>(-.8,.6);
    const TinyVector<double,2> reflection(1.,0.);
    const auto energy = [&](const auto &a, const auto &b) {
        // Literal independent H=-x^T J y, with the original ordered endpoints.
        return -(a.state_[0]*b.state_[0] + .7*a.state_[0]*b.state_[1]
                 -.4*a.state_[1]*b.state_[0] -.2*a.state_[1]*b.state_[1]);
    };
    for (const auto vertex : {left,right}) {
        const auto incident = *boost::out_edges(vertex,graph.graph()).first;
        const auto oriented = spinmc_local_bond_matrix(matrix,graph.graph(),incident,sources);
        auto next_x=x,next_y=y;
        if (vertex==left) next_x.update(reflection); else next_y.update(reflection);
        const double expected = energy(next_x,next_y)-energy(x,y);
        const double actual = vertex==left ? bond_energy_change(x,y,oriented,reflection)
                                          : bond_energy_change(y,x,oriented,reflection);
        EXPECT_NEAR(actual,expected,1.e-15);
    }
}

TEST(SpinMCScientificCorrections, OrientedXYExchangeMatchesIndependentRotatedDimerPartition) {
    auto input=parameters("XY");
    input["S"]=1.; input["THERMALIZATION"]=1000;
    std::ostringstream matrix;
    // J=.6R(.4). Rotating one integration variable leaves the exact XY
    // partition Z proportional to I0(beta*.6), including nonsymmetric J.
    const double diagonal=.6*std::cos(.4), off_diagonal=.6*std::sin(.4);
    matrix<<std::setprecision(17)<<diagonal<<' '<<-off_diagonal<<' '
          <<off_diagonal<<' '<<diagonal;
    input["J"]=matrix.str();
    original_run run(input);
    for (unsigned warmup=0; warmup<1000; ++warmup) run.worker->run();
    for(unsigned update=0;update<100000;++update) {
        run.worker->run();
    }
    const double correlation=boost::math::cyl_bessel_i(1,.5)/boost::math::cyl_bessel_i(0,.5);
    EXPECT_NEAR(alps::RealObsevaluator(run.measurements->get_measurements()["Energy"]).mean(),-.6*correlation,.025);
    EXPECT_NEAR(alps::RealObsevaluator(run.measurements->get_measurements()["Energy^2"]).mean(),.36*(1-2*correlation),.025);
}

TEST(SpinMCScientificCorrections, SelfBondUsesOneSymmetricQuadraticReflectionDifference) {
    const auto input=parameters("XY");
    std::string text="1 .7 -.4 -.2";
    const auto symmetric=spinmc_self_bond_matrix(SquareMatrix<double,2>(text,input));
    ONMoment<2> spin; spin.state_=TinyVector<double,2>(.6,.8);
    const TinyVector<double,2> reflection(1.,0.);
    // H=-x²-.3xy+.2y². The antisymmetric exchange drops out on a self-bond.
    EXPECT_NEAR(onsite_energy_change(spin,symmetric,reflection),.288,1.e-15);
    MIdMatrix<double,2> isotropic(.6);
    EXPECT_DOUBLE_EQ(onsite_energy_change(spin,isotropic,reflection),0.);
}
