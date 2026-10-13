// SPDX-License-Identifier: MIT
// Unit tests for spinmc::resolve_couplings. Terms resolve from parameters and
// the site types and bond classes of a lattice, without building the lattice.
#include "couplings.h"
// matrices.h relies on alps::evaluate being declared before it.
#include "matrices.h"

#include <alps/expression.h>
#include <gtest/gtest.h>

#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using spinmc::bond_class;
using spinmc::resolve_couplings;

alps::Parameters parameters(const std::string& text) {
    alps::Parameters p;
    std::istringstream input(text);
    p.parse(input);
    return p;
}

const std::vector<bond_class> no_bonds;

TEST(SpinmcCouplings, FallsBackFromTypedToGenericToDefaultSpin) {
    const auto c = resolve_couplings<MIdMatrix<double, 3>>(parameters("S1=2.5;S=1.5"), {0, 1},
                                                           no_bonds, false, "1.0");
    ASSERT_EQ(c.sites.size(), 2u);
    EXPECT_EQ(c.sites[0].type, 0);
    EXPECT_EQ(c.sites[0].spin_parameter, "S");
    EXPECT_DOUBLE_EQ(c.sites[0].spin, 1.5);
    EXPECT_DOUBLE_EQ(c.sites[0].spinfactor, 1.5);
    EXPECT_EQ(c.sites[1].spin_parameter, "S1");
    EXPECT_DOUBLE_EQ(c.sites[1].spinfactor, 2.5);

    const auto classical = resolve_couplings<MIdMatrix<double, 3>>(parameters("J=1"), {3},
                                                                   no_bonds, false, "1.0");
    EXPECT_EQ(classical.sites[0].spin_parameter, "");
    EXPECT_DOUBLE_EQ(classical.sites[0].spin, 1.);
    EXPECT_DOUBLE_EQ(classical.sites[0].spinfactor, 1.);

    const auto quantum = resolve_couplings<MIdMatrix<double, 3>>(parameters("J=1"), {3},
                                                                 no_bonds, true, "1.0");
    EXPECT_DOUBLE_EQ(quantum.sites[0].spin, 0.5);
    EXPECT_DOUBLE_EQ(quantum.sites[0].spinfactor, std::sqrt(0.75));
}

TEST(SpinmcCouplings, ScalesOnsiteTermsBySpinSquared) {
    const auto c = resolve_couplings<DiagMatrix<double, 3>>(
        parameters(R"p(D1="0.1 0.2 0.3";S=2)p"), {0, 1}, no_bonds, false, "1.0");
    EXPECT_TRUE(c.general_case);
    EXPECT_EQ(c.sites[0].onsite_parameter, "");
    EXPECT_EQ(c.sites[1].onsite_parameter, "D1");
    for (int i = 0; i < 3; ++i) {
        EXPECT_DOUBLE_EQ(c.sites[0].selfinteraction.getElement(i, i), 0.);
        EXPECT_DOUBLE_EQ(c.sites[1].selfinteraction.getElement(i, i), 0.4 * (i + 1));
    }
    EXPECT_FALSE((resolve_couplings<DiagMatrix<double, 3>>(parameters("S=2"), {0, 1}, no_bonds,
                                                           false, "1.0")
                      .general_case));
}

TEST(SpinmcCouplings, RejectsNonSymmetricOnsiteTerms) {
    try {
        resolve_couplings<SquareMatrix<double, 3>>(parameters(R"p(D="1 2 3 4 5 6 7 8 9")p"),
                                                   {0}, no_bonds, false, "1.0");
        FAIL() << "expected an exception";
    } catch (const std::runtime_error& error) {
        EXPECT_STREQ(error.what(), "Invalid matrix for D");
    }
}

// J<type> wins over the default text, which is the value of J or "1.0".
// Couplings change sign in the classical convention.
TEST(SpinmcCouplings, ResolvesCouplingsPerBondClass) {
    const std::vector<bond_class> classes = {{0, 0, 0}, {1, 0, 0}, {0, 0, 1}};
    const auto c = resolve_couplings<MIdMatrix<double, 1>>(parameters("J0=2;J=0.5"), {0, 1},
                                                           classes, false, "0.5");
    ASSERT_EQ(c.bonds.size(), 3u);
    EXPECT_EQ(c.bonds[0].coupling_parameter, "J0");
    EXPECT_DOUBLE_EQ(c.bonds[0].coupling.getElement(0, 0), -2.);
    EXPECT_EQ(c.bonds[1].bonds.type, 1);
    EXPECT_EQ(c.bonds[1].coupling_parameter, "J");
    EXPECT_DOUBLE_EQ(c.bonds[1].coupling.getElement(0, 0), -0.5);
    EXPECT_EQ(c.bonds[2].bonds.target_type, 1);
    EXPECT_EQ(c.bonds[2].coupling_parameter, "J0");

    const auto unset = resolve_couplings<MIdMatrix<double, 1>>(parameters("S=1"), {0},
                                                               {{0, 0, 0}}, true, "1.0");
    EXPECT_EQ(unset.bonds[0].coupling_parameter, "");
    EXPECT_DOUBLE_EQ(unset.bonds[0].coupling.getElement(0, 0), 1.);
}

// The model library names the couplings of bond types 1 and 2 J' and J''.
// They are not read, so bonds of those types without J1 or J2 are rejected
// instead of taking J. J' stays usable in expressions.
TEST(SpinmcCouplings, RejectsPrimedCouplingsThatWouldBeIgnored) {
    const std::vector<bond_class> legs = {{0, 0, 0}};
    const std::vector<bond_class> ladder = {{0, 0, 0}, {1, 0, 0}};
    const std::vector<bond_class> second = {{0, 0, 0}, {2, 0, 0}};
    const auto error = [](const std::string& text, const std::vector<bond_class>& classes) {
        try {
            resolve_couplings<MIdMatrix<double, 1>>(parameters(text), {0}, classes, false, "1");
        } catch (const std::runtime_error& e) {
            return std::string(e.what());
        }
        return std::string();
    };
    const std::string ignored_j1 =
        R"(J' is not read: bond type 1 takes J1, falling back to J. Set J1, for example J1="J'", or remove J')";
    EXPECT_EQ(error("J=1;J'=0.5", ladder), ignored_j1);
    EXPECT_EQ(error("J'=0.5", ladder), ignored_j1);
    EXPECT_EQ(error("J=1;J''=0.5", second),
              R"(J'' is not read: bond type 2 takes J2, falling back to J. Set J2, for example J2="J''", or remove J'')");

    // Without bonds of the primed type, or with their J<type>, J' is accepted.
    EXPECT_EQ(error("J=1;J'=0.5", legs), "");
    EXPECT_EQ(error("J=1;J''=0.5", ladder), "");
    EXPECT_EQ(error("J=1;J'=0.5;J1=0.25", ladder), "");
    const auto c = resolve_couplings<MIdMatrix<double, 1>>(parameters(R"p(J=1;J'=0.25;J1="J'")p"),
                                                           {0}, ladder, false, "1");
    EXPECT_EQ(c.bonds[1].coupling_parameter, "J1");
    EXPECT_DOUBLE_EQ(c.bonds[0].coupling.getElement(0, 0), -1.);
    EXPECT_DOUBLE_EQ(c.bonds[1].coupling.getElement(0, 0), -0.25);
}

TEST(SpinmcCouplings, ClassifiesCouplingSigns) {
    const std::vector<bond_class> one = {{0, 0, 0}};
    const auto antiferro = resolve_couplings<MIdMatrix<double, 1>>(parameters("J=1"), {0}, one,
                                                                   false, "1");
    EXPECT_LT(antiferro.bonds[0].det, 0.);
    EXPECT_TRUE(antiferro.ferromagnetic);
    EXPECT_FALSE(antiferro.antiferromagnetic);
    const auto ferro = resolve_couplings<MIdMatrix<double, 1>>(parameters("J=-1"), {0}, one,
                                                               false, "-1");
    EXPECT_FALSE(ferro.ferromagnetic);
    EXPECT_TRUE(ferro.antiferromagnetic);
}

TEST(SpinmcCouplings, LeavesParametersUnchanged) {
    const auto p = parameters("J=1");
    resolve_couplings<MIdMatrix<double, 3>>(p, {0, 1}, {{0, 0, 1}}, false, "1");
    EXPECT_EQ(p.size(), 1u);
}

// Current behavior, pinned so that a change is deliberate: terms are drawn
// from the disorder generator site type by site type, S before D, and then
// bond class by bond class.
TEST(SpinmcCouplings, DrawsRandomTermsInResolutionOrder) {
    alps::Disorder::seed(5);
    std::vector<double> draws;
    for (int i = 0; i < 5; ++i)
        draws.push_back(alps::Disorder::random());

    alps::Disorder::seed(5);
    const auto c = resolve_couplings<MIdMatrix<double, 1>>(
        parameters(R"p(S="random()";D="random()";J1="random()")p"), {0, 1}, {{1, 0, 1}}, false,
        "1.0");
    EXPECT_DOUBLE_EQ(c.sites[0].spin, draws[0]);
    EXPECT_DOUBLE_EQ(c.sites[0].selfinteraction.getElement(0, 0),
                     draws[1] * draws[0] * draws[0]);
    EXPECT_DOUBLE_EQ(c.sites[1].spin, draws[2]);
    EXPECT_DOUBLE_EQ(c.sites[1].selfinteraction.getElement(0, 0),
                     draws[3] * draws[2] * draws[2]);
    EXPECT_DOUBLE_EQ(c.bonds[0].coupling.getElement(0, 0), -draws[4]);
}

} // namespace
